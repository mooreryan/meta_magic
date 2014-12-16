#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

begin
  require 'trollop'
  require 'shell/executer.rb'
  require 'parse_fasta'
  require 'fileutils'
rescue LoadError => e
  bad_file = e.message.sub(/^cannot load such file -- /, '')
  abort("ERROR: #{e.message}\nTry running: gem install #{bad_file}")
end

opts = Trollop.options do
  banner <<-EOS

  Steps:

  1) Count reads
  2) Interleave reads
  3) Get quality stats on starting sequences
  4) Flash reads together
  5) Quality filter
  6) Separate into pe and se files
  7) Count reads in pe and se files
  8) Get quality stats on final sequences

  Options:
  EOS
  opt(:left, 'Left reads (.fq.gz)', type: :string)
  opt(:right, 'Right reads (.fq.gz)', type: :string)
  opt(:prefix, 'Prefix for output', type: :string)
  opt(:quality, 'Minimum quality score to keep',
      type: :int,
      default: 30)
  opt(:percent,
      'Minimum percent of bases that must have [-q] quality',
      type: :int,
      default: 50)
  opt(:outdir, 'Output directory', type: :string)
  opt(:threads,
      'Number of threads to use ' +
      '(if greater than one, will use pigz)',
      type: :int,
      default: 1)

  opt(:kmer, 'Diginorm kmer size', type: :int, default: 20)
  opt(:coverage, 'Diginorm max cov', type: :int, default: 10)
  opt(:ram, 'How much ram to use for diginorm hash',
      type: :int, default: 4)
end

if opts[:left].nil?
  Trollop.die :left, "You must enter a file name"
elsif !File.exist? opts[:left]
  Trollop.die :left, "The file must exist"
end

if opts[:right].nil?
  Trollop.die :right, "You must enter a file name"
elsif !File.exist? opts[:right]
  Trollop.die :right, "The file must exist"
end

if opts[:outdir].nil?
  Trollop.die :outdir, "You must enter an output directory"
elsif !File.exist? opts[:outdir]
  begin
    FileUtils.mkdir(opts[:outdir])
  rescue EACCES => e
    abort("ERROR: #{e.message}\n\n" +
          "It appears you don't have the proper permissions")
  end
end

if opts[:prefix].nil?
  Trollop.die :prefix, "Don't forget to specify a prefix!"
end

if opts[:threads] < 1
  Trollop.die :threads, "Threads must be positive"
end

if opts[:ram] < 4
  Trollop.die :ram, "You must specify at least 4 gb of ram"
end

ram = "#{opts[:ram] / 4.0}e9"

######################################################################
##### functions ######################################################
######################################################################

def run_it(cmd)
  begin
    $stderr.puts "\nRUNNING: #{cmd}"
    cmd_outerr = Shell.execute!(cmd)
    $stdout.puts cmd_outerr.stdout unless cmd_outerr.stdout.empty?
    $stderr.puts cmd_outerr.stderr unless cmd_outerr.stderr.empty?
    
    return cmd_outerr
  rescue RuntimeError => e
    # print stderr if bad exit status
    abort(e.message)
  end
end

def parse_fname(fname)
  { dir: File.dirname(fname),
    base: File.basename(fname, File.extname(fname)),
    ext: File.extname(fname) }
end

def elapsed_time(old_time)
  (Time.now - old_time).quo(60).round(2)
end

def rand_str(len=10)
  len.times.map { ('a'..'z').to_a.sample }.join
end

# quality info before filtering
def qual_stats(interleaved_reads, starting=false)
  home = '/home/moorer'
  fastx = 'vendor/fastx_toolkit-0.0.14/bin'

  if starting
    output = "#{interleaved_reads}.pre_qc.stats.txt"
  else
    output = "#{interleaved_reads}.stats.txt"
  end
  
  cmd =
    "#{home}/#{fastx}/fastx_quality_stats -i #{interleaved_reads} " +
    "-o #{output}"
  run_it(cmd)
end

countfq = lambda do |infile, threads|
  "cat #{infile} | echo $((`wc -l`/4))"
end

zcountfq = lambda do |infile, threads|
  if threads == 1
    "zcat #{infile} | echo $((`wc -l`/4))"
  elsif threads > 1
    "unpigz -d -c -p #{threads} #{infile} | echo $((`wc -l`/4))"
  end
end

def count_reads(*files, out_fname: '', count_fn: '', threads: 1)
  counts = files.map do |file|
    run_it(count_fn.call(file, threads)).stdout.to_i
  end
  
  File.open(out_fname, 'w') do |f|
    counts.each_with_index do |count, idx|
      f.puts [files[idx], count].join(' ')
    end
  end
end

######################################################################
##### pipeline #######################################################
######################################################################

# source proper env vars
run_it(". /home/moorer/.bash_profile")


#### scripts #########################################################

home = '/home/moorer'
khmer = 'vendor/khmer/scripts'
sandbox = 'vendor/khmer/sandbox'
fastx = 'vendor/fastx_toolkit-0.0.14/bin'
interleave = "#{home}/#{khmer}/interleave-reads.py"
q_filter = "#{home}/#{fastx}/fastq_quality_filter"
extract_paired_reads = "#{home}/#{khmer}/extract-paired-reads.py"
flash = '/usr/local/bin/flash'
readstats = "#{home}/vendor/khmer/sandbox/readstats.py"

diginorm = {
  normalize_by_median: "#{home}/#{khmer}/normalize-by-median.py",
  filter_below_abund: "#{home}/#{sandbox}/my-filter-below-abund.py",
}

# make sure they are there
should_exit = false
programs = [interleave,
            q_filter,
            extract_paired_reads,
            flash,
            readstats]
programs.each do |program|
  unless File.exist?(program)
    $stderr.puts "ERROR: #{program} doesn't exist!"
    should_exit = true
  end
end
exit if should_exit

diginorm.values.each do |program|
  unless File.exist?(program)
    $stderr.puts "ERROR: #{program} doesn't exist!"
    should_exit = true
  end
end
exit if should_exit

#### file names ######################################################

interleaved = "#{opts[:prefix]}.interleaved"

original_read_counts_fname =
  File.join(opts[:outdir], "#{opts[:prefix]}.original.counts.txt")

interleaved_reads_fname =
  File.join(opts[:outdir], "#{interleaved}.fq")

old_flashed_fname =
  File.join(opts[:outdir], "flashed.extendedFrags.fastq")
flashed_fname =
  File.join(opts[:outdir], "#{interleaved}.flashed.fq")
old_nonflashed_fname =
  File.join(opts[:outdir], "flashed.notCombined.fastq")
nonflashed_fname =
  File.join(opts[:outdir], "#{interleaved}.nonflashed.fq")

# PREFIX.interleaved.flashed.filtered.fq
flashed_filtered_reads_fname =
  File.join(opts[:outdir],
            "#{parse_fname(flashed_fname)[:base]}.filtered.fq")
# PREFIX.interleaved.nonflashed.filtered.fq
nonflashed_filtered_reads_fname =
  File.join(opts[:outdir],
            "#{parse_fname(nonflashed_fname)[:base]}.filtered.fq")

# contains all se reads: flashed and qc orphans
# PREFIX.flashed_and_filtered.se.fq
flashed_and_filtered_se_fname =
  File.join(opts[:outdir],
            "#{opts[:prefix]}.flashed_and_filtered.se.fq")
nonflashed_and_filtered_pe_fname =
  "#{nonflashed_filtered_reads_fname}.pe"

pe_se_counts_fname =
  File.join(opts[:outdir], "#{opts[:prefix]}.pe_se.counts.txt")

pe_fname =
  File.join(opts[:outdir], "#{opts[:prefix]}.pe.filtered.fq")
se_fname =
  File.join(opts[:outdir], "#{opts[:prefix]}.se.filtered.fq")

pe_gz_fname =
  File.join(opts[:outdir], "#{opts[:prefix]}.pe.filtered.fq.gz")
se_gz_fname =
  File.join(opts[:outdir], "#{opts[:prefix]}.se.filtered.fq.gz")

info_dir = File.join(opts[:outdir], 'info')
stats_dir = File.join(info_dir, 'stats')
counts_dir = File.join(info_dir, 'counts')
flash_dir = File.join(info_dir, 'flash_info')

# diginorm stuff
digi_dir = File.join(opts[:outdir], 'digi_dir')
FileUtils.mkdir(digi_dir)
hash_table_fname =
  File.join(digi_dir,
            "#{opts[:prefix]}.c#{opts[:coverage]}.k#{opts[:kmer]}.kh")
pe_keep_fname =
  File.join(digi_dir,
            "#{parse_fname(pe_gz_fname)[:base]}.gz.keep")
se_keep_fname =
  File.join(digi_dir,
            "#{parse_fname(se_gz_fname)[:base]}.gz.keep")

#### count reads in each file ########################################

cmd =
  "#{readstats} #{opts[:left]} #{opts[:right]}" +
  "> #{original_read_counts_fname}"
run_it(cmd)

#### interleave reads ################################################

cmd =
  "#{interleave} -o #{interleaved_reads_fname} #{opts[:left]} " +
  "#{opts[:right]}"
run_it(cmd)

#### quality stats ###################################################

qual_stats(interleaved_reads_fname, true)

#### flash ###########################################################

# flash reads
cmd =
  "#{flash} --interleaved --output-prefix flashed " +
  "--output-directory #{opts[:outdir]} --threads #{opts[:threads]} " +
  "#{interleaved_reads_fname}"
run_it(cmd)

# fix naming scheme
FileUtils.mv(old_flashed_fname, flashed_fname)
FileUtils.mv(old_nonflashed_fname, nonflashed_fname)

# remove old interleaved file
FileUtils.rm(interleaved_reads_fname)

#### quality filter ##################################################

# do the flashed se reads
cmd =
  "#{q_filter} -Q33 -q #{opts[:quality]} -p #{opts[:percent]} " +
  "-i #{flashed_fname} > #{flashed_filtered_reads_fname}"
run_it(cmd)

cmd =
  "#{q_filter} -Q33 -q #{opts[:quality]} -p #{opts[:percent]} " +
  "-i #{nonflashed_fname} > #{nonflashed_filtered_reads_fname}"
run_it(cmd)

FileUtils.rm(flashed_fname)
FileUtils.rm(nonflashed_fname)

#### separate properly paired and orfaned reads ######################

# will output:
#   #{interleaved}.filtered.fq.pe <-- proper pairs
#   #{interleaved}.filtered.fq.se <-- orphaned reads

cmd = "#{extract_paired_reads} #{nonflashed_filtered_reads_fname}"
run_it(cmd)

# move extracted files to output folder
FileUtils.mv(Dir.glob("#{interleaved}*"), opts[:outdir])
FileUtils.rm(nonflashed_filtered_reads_fname)

# combine the flashed reads (pe) and the qc-ed pe reads into one pe
# file
cmd =
  "cat #{flashed_filtered_reads_fname} " +
  "#{nonflashed_filtered_reads_fname}.se " +
  "> #{flashed_and_filtered_se_fname}"
run_it(cmd)
FileUtils.rm(flashed_filtered_reads_fname)
FileUtils.rm("#{nonflashed_filtered_reads_fname}.se")

#### quality stats ###################################################

# use a sane name
FileUtils.mv(nonflashed_and_filtered_pe_fname, pe_fname)
FileUtils.mv(flashed_and_filtered_se_fname, se_fname)

qual_stats(pe_fname)
qual_stats(se_fname)

#### gzip the pe and se files ########################################

if opts[:threads] == 1
  cmd = "gzip #{pe_fname}"
  run_it(cmd)

  cmd = "gzip #{se_fname}"
  run_it(cmd)
else
  cmd =
    "pigz --best -p #{opts[:threads]} #{pe_fname}"
  run_it(cmd)

  cmd =
    "pigz --best -p #{opts[:threads]} #{se_fname}"
  run_it(cmd)
end

#### count sequences in pe and se files ##############################

cmd =
  "#{readstats} #{pe_gz_fname} #{se_gz_fname} > #{pe_se_counts_fname}"
run_it(cmd)

#### move extra output to its own folder #############################

FileUtils.mkdir([info_dir, stats_dir, counts_dir, flash_dir])
FileUtils.mv(Dir.glob(File.join(opts[:outdir], '*hist*')), flash_dir)

FileUtils.mv(Dir.glob(File.join(opts[:outdir], '*stats*txt')),
             stats_dir)
FileUtils.mv(Dir.glob(File.join(opts[:outdir], '*counts*txt')),
             counts_dir)


#### diginorm ########################################################

# the paired end reads first
cmd =
  "#{diginorm[:normalize_by_median]} -k #{opts[:kmer]} " +
  "-C #{opts[:coverage]} -N 4 -x #{ram} " +
  "-o #{pe_keep_fname} " +
  "--savetable #{hash_table_fname} -p #{pe_gz_fname}"
run_it(cmd)

# now the se reads
cmd =
  "#{diginorm[:normalize_by_median]} " +
  "-C #{opts[:coverage]} " +
  "-o #{se_keep_fname} --loadtable #{hash_table_fname} " +
  "--savetable #{hash_table_fname} #{se_gz_fname}"
run_it(cmd)

#### filter below abund ############################################

cmd =
  "#{diginorm[:filter_below_abund]} " +
  "#{opts[:threads]} " +
  "#{opts[:coverage] * 5} " +
  "#{hash_table_fname} " +
  "#{pe_keep_fname} " +
  "#{se_keep_fname}"
run_it(cmd)

FileUtils.mv(Dir.glob("*keep.below"), digi_dir)

puts
puts "Done!"
