#!/usr/bin/env python3

import csv
import os

#############
# FUNCTIONS #
#############

# lookup wildcards for merge function
def merge_wildcard_resolver(wildcards):
    '''
    The merge_per_sample rule passes two wildcards, `sample_name` and `read`.
    Look up the file names using the name_to_sample dict and return a list of
    files that match by filename and read number
    '''
    my_filename = name_to_sample[wildcards.sample_name]
    return sorted(set(x for x in fastq_files
                  if (my_filename in x
                      and wildcards.read in x)))


###########
# GLOBALS #
###########

read_dir = 'data/reads'
sample_key = 'data/NZGL01923A_C9KB3ANXX_sequencing_summary.txt'
all_samples = ['4872EtOH-1', '4872QMP-2', '48EtOH-2', '4872EtOH-2',
               '4872QMP-1', '12EtOH-2', '48QMP-1', '24QMP-1',
               '48EtOH-1', '12EtOH-1', '4872stay-2', '4872stay-1',
               '48QMP-2', '12QMP-2', '24QMP-2', '24EtOH-2', '12QMP-1',
               '24EtOH-1']


# find fastq files
read_dir_files = [(dirpath, filenames) for (dirpath, dirnames, filenames)
                  in os.walk(read_dir)]
fastq_files = []
for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            fastq_files.append(os.path.join(dirpath, filename))

# generate name to filename dictionary
with open(sample_key) as csvfile:
    csvreader = csv.reader(csvfile, delimiter="\t")
    next(csvreader)
    name_to_sample = {x[1]: x[4] for x in csvreader
                      if x[1] in all_samples}


#########
# RULES #
#########

# TODO: write logs and stats to a proper log directory

rule all:
    input:
        expand('output/bbduk/{sample_name}_{read}.fastq.gz',
               sample_name=all_samples, read=['R1', 'R2'])

rule merge_per_sample:
    threads:
        1
    input:
        merge_wildcard_resolver
    output:
        'output/merged/{sample_name}_{read}.fastq.gz'
    shell:
        'cat {input} '
        '> {output}'

rule trim_decon:
    threads:
        18
    input:
        r1 = 'output/merged/{sample_name}_R1.fastq.gz',
        r2 = 'output/merged/{sample_name}_R2.fastq.gz',
    output:
        r1 = 'output/bbduk/{sample_name}_R1.fastq.gz',
        r2 = 'output/bbduk/{sample_name}_R2.fastq.gz',
        stats = 'output/bbduk/{sample_name}_stats.txt'
    shell:
        'bin/bbmap/bbduk.sh '
        'threads={threads} '
        '-Xmx100g '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'ref=bin/bbmap/resources/adapters.fa ' 
        '| bin/bbmap/bbduk.sh '
        'threads={threads} '
        '-Xmx100g '
        'in=stdin.fastq '
        'out={output.r1} '
        'out2={output.r2} '
        'ref=bin/bbmap/resources/sequencing_artifacts.fa.gz '
        'k=31 hdist=1 stats={output.stats}'
