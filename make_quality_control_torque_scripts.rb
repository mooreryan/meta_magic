#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

require 'fileutils'

torque_dir = File.absolute_path(ARGV.first)
FileUtils.mkdir(torque_dir) unless File.exist?(torque_dir)

FileUtils.cd('/data/moorer/repository/cow/viral_metagenomes')

Dir.glob('Sample*').each do |folder|
  left, right = Dir.glob("#{folder}/*gz").map { |f| File.absolute_path(f) }
  prefix = folder
  outdir = File.join(File.absolute_path(folder), 'qc')
  threads = 8

  script_name = File.join(torque_dir, "#{prefix}.quality_control.submitter.qs")
  puts script_name
  File.open(script_name, 'w') do |f|
    f.puts(
"#!/bin/bash

#PBS -N Jurassic5-QC
#PBS -l walltime=4:00:00,nodes=biohen36:ppn=#{threads},cput=32:00:00
#PBS -d /home/moorer/runt
#PBS -e /home/moorer/runt/oe
#PBS -o /home/moorer/runt/oe

## code to run here

hostname
date

script=/home/moorer/git_repos/pipelines/metagenomics/process_illumina_metagenome/quality_control.rb

left=#{left}
right=#{right}
prefix=#{prefix}
outdir=#{outdir}
threads=#{threads}

time ruby $script -l $left -r $right -p $prefix -o $outdir -t $threads

date
echo 'done!'
")
  end
end
