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

  Options:
  EOS
  opt(:left, 'Left reads (.fq.gz)', type: :string)
  opt(:right, 'Right reads (.fq.gz)', type: :string)
  opt(:prefix, 'Prefix for output', type: :string)
  opt(:quality, 'Minimum quality score to keep', type: :int, default: 30)
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
  opt(:print_only, 'Print commands w/o calling them', type: :boolean,
      default: false)

  # opt(:bin, 'The directory contianing all the scripts',
  #     default: '/Users/moorer/projects/cone_janws/scripts')
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
  Trollop.die :outdir, "The directory must exist"
end

if opts[:prefix].nil?
  Trollop.die :prefix, "Don't forget to specify a prefix!"
end

if opts[:threads] < 1
  Trollop.die :threads, "Threads must be positive"
end

  
# if !File.exist? opts[:bin]
#   Trollop.die :bin, "The file must exist"
# end

######################################################################
##### functions ######################################################
######################################################################

def run_it(cmd, print_only)
  if print_only
    $stderr.puts cmd
  else
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

def zcountfq(infile)
  "zcat #{infile} | echo $((`wc -l`/4))"
end

# quality info before filtering
def qual_stats(interleaved_reads, print_only)
  home = '/home/moorer'
  fastx = 'vendor/fastx_toolkit-0.0.14/bin'
  
  cmd = "#{home}/#{fastx}/fastx_quality_stats -i #{interleaved_reads} -o #{interleaved_reads}.stats.txt"
  run_it(cmd, print_only)

  # cmd = "#{home}/#{fastx}/fastq_quality_boxplot_graph.sh -i #{interleaved_reads}.stats.txt -o #{interleaved_reads}.quality.png -t #{prefix}"
  # run_it(cmd, opts[:print_only])

  # cmd = "#{home}/#{fastx}/fastx_nucleotide_distribution_graph.sh -i #{interleaved_reads}.stats.txt -o #{interleaved_reads}.nuc_dist.png -t #{prefix}"
  # run_it(cmd, opts[:print_only])
end

######################################################################
##### pipeline #######################################################
######################################################################

# source proper env vars
run_it(". /home/moorer/.bash_profile", opts[:print_only])

home = '/home/moorer'
khmer = 'vendor/khmer/scripts'
fastx = 'vendor/fastx_toolkit-0.0.14/bin'
zcountfq = 
interleave = "#{home}/#{khmer}/interleave-reads.py"
q_filter = "#{home}/#{fastx}/fastq_quality_filter"
extract_paired_reads = "#{home}/#{khmer}/extract-paired-reads.py"

combined = "#{opts[:prefix]}_combined"

# count reads in each file

if opts[:print_only]
  num_sequences_left = 10
  num_sequences_right = 10
else
  num_sequences_left = run_it(zcountfq(opts[:left]),
                              opts[:print_only]).stdout.to_i
  num_sequences_right = run_it(zcountfq(opts[:right]),
                               opts[:print_only]).stdout.to_i
end  
counts = "#{combined}.counts.txt"
File.open(File.join(opts[:outdir], counts), 'w') do |f|
  f.puts [opts[:left], num_sequences_left].join(' ')
  f.puts [opts[:right], num_sequences_right].join(' ')
end

# interleave reads
interleaved_reads = File.join(opts[:outdir], "#{combined}.fq")
cmd = "#{interleave} -o #{interleaved_reads} #{opts[:left]} #{opts[:right]}"
run_it(cmd, opts[:print_only])

qual_stats(interleaved_reads, opts[:print_only])

# qualty filter
filtered_reads = File.join(opts[:outdir], "#{combined}.filtered.fq") 
cmd = "#{q_filter} -Q33 -q #{opts[:quality]} -p #{opts[:percent]} -i #{interleaved_reads} > #{filtered_reads}"
run_it(cmd, opts[:print_only])

if opts[:print_only]
  $stderr.puts "FileUtils.rm(interleaved_reads) if File.exist?(interleaved_reads)"
else
  FileUtils.rm(interleaved_reads) if File.exist?(interleaved_reads)
end

qual_stats(filtered_reads, opts[:print_only])

# separate properly paired and orfaned reads
# will output:
#   #{combined}.filtered.fq.pe <-- proper pairs
#   #{combined}.filtered.fq.se <-- orphaned reads
cmd = "#{extract_paired_reads} #{filtered_reads}"
run_it(cmd, opts[:print_only])
# moved extracted files to output folder

if opts[:print_only]
  $stderr.puts("FileUtils.mv(Dir.glob(\"#{combined}*\"), opts[:outdir])\n" +
               "FileUtils.rm(filtered_reads) if File.exist?(filtered_reads)")
else
  FileUtils.mv(Dir.glob("#{combined}*"), opts[:outdir])
  FileUtils.rm(filtered_reads) if File.exist?(filtered_reads)
end
# gzip the pe and se files
pe_out = File.join(opts[:outdir], "#{opts[:prefix]}.pe.filtered.fq.gz")
se_out = File.join(opts[:outdir], "#{opts[:prefix]}.se.filtered.fq.gz")
if opts[:threads] == 1
  cmd = "gzip -c #{filtered_reads}.pe > #{pe_out}"
  run_it(cmd, opts[:print_only])

  cmd = "gzip -c #{filtered_reads}.se > #{se_out}"
  run_it(cmd, opts[:print_only])
else
  cmd =
    "pigz -c --best -p #{opts[:threads]} #{filtered_reads}.pe > #{pe_out}"
  run_it(cmd, opts[:print_only])

  cmd =
    "pigz -c --best -p #{opts[:threads]} #{filtered_reads}.se > #{se_out}"
  run_it(cmd, opts[:print_only])
end

if opts[:print_only]
  $stderr.puts "FileUtils.rm(Dir.glob(\"#{filtered_reads}.?e\"))"
else
  FileUtils.rm(Dir.glob("#{filtered_reads}.?e"))
end

# move extra output to its own folder
info_dir = File.join(opts[:outdir], 'info')
unless File.exist?(info_dir)

  if opts[:print_only]
    $stderr.puts "info_dir = FileUtils.mkdir(info_dir)"
  else
    info_dir = FileUtils.mkdir(info_dir)
  end
end

if opts[:print_only]
  $stderr.puts("FileUtils.mv(Dir.glob(File.join(opts[:outdir], '*.txt')), info_dir)")
else
  FileUtils.mv(Dir.glob(File.join(opts[:outdir], '*.txt')), info_dir)
end

puts
