'''

    Utilities for analyzing ribosome profiling data, storing ribosome density files
    
    Copyright (C) 2021  Allen Buskirk, buskirk@jhmi.edu
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
'''


from BCBio import GFF
import csv
from Bio import Seq
import struct
from multiprocessing import Process
import os
import cPickle as pickle
from datetime import datetime



def multiprocess(function, arguments, threads):
    ''' 
    Multiprocess utility to run multiple functions in parallel
    Useful for creating density or running analysis on multiple files
    '''
    processes = []
    
    for arg in arguments:
        p = Process(target = function, args = arg)
        processes.append(p)
        p.start()
    
    for p in processes:
        p.join()
        
    return      
  
def countstowig(countsfile,filestring):
    import random
    f=open(filestring+".wig","w")
    filestring=filestring.partition("_")[0][-3:]

    random.seed(filestring)
    c1=random.randint(0,255)
    random.seed(filestring+"xxx")
    c2=random.randint(0,255)
    random.seed(filestring+"000")
    c3=random.randint(0,255)

    f.write("track type=wiggle_0 name=tracklabel viewLimits=-5:5 color="+str(c1)+','+str(c2)+','+str(c3)+"\n")
    for chrom in countsfile.keys():
        if chrom[0:3]=='chr':
            f.write("fixedStep  chrom="+chrom+"  start=1  step=1\n")
        else:
            f.write("fixedStep  chrom=\""+chrom+"\"  start=1  step=1\n")

        for position in countsfile[chrom]:
            f.write(str(position)+"\n")
    f.close()
    
def makePickle(data, path_pickle, protocol=pickle.HIGHEST_PROTOCOL):
    f = open(path_pickle, 'w')
    pickle.dump(data, f, protocol=protocol)
    f.close()
    
    
def unPickle(path_pickle):
    #returns the pickled object stored in a pickle file
    f = open(path_pickle, 'r')
    data = pickle.load(f)
    f.close()
    return data


def density_3(fname, chr_sam, minlength, maxlength, path_wig, path_den, path_gff, source):
    
    fname = fname
    chr_sam = chr_sam
    minlength = minlength
    maxlength = maxlength
    GFFgen = GFF.parse(path_gff)

    # open chr aligned sam file
    f_samfile = open(chr_sam)
    samfile = csv.reader(f_samfile,delimiter='	')
    
    # dictionaries to hold read counts
    density_plus = {}
    density_minus = {}
    
    if minlength < 0 or maxlength < 0:
        print "Error. Length input not valid."
        return(0)
    
    # Makes 2 sets of indices, one for all reads, and another for size separated:
    for sequence in GFFgen:
        density_plus[sequence.id]  = [0 for x in range(len(sequence))]
        density_minus[sequence.id] = [0 for x in range(len(sequence))]
   
    total_reads = 0
    mapped_reads = 0

    # Loop through the samfile.
    for read in samfile:
        if read[0][0] == '@':   # Ignore header lines.
            continue

        if read[1] == '4':      # A bowtie mismatch.  
            continue

        chrom = read[2]             # chromosome identified for read in bowtie
        readid = read[0]            # read id
        startp = int(read[3]) -1    # start position. Need to subtract 1 since genomic sequence starts at 1,
        seq = Seq.Seq(read[9])      # sequence of the read
        length = len(seq)           # length of read
        
        if chrom not in density_plus.keys():
            print "Error: Bowtie index and GFF do not match"
        
        total_reads += 1

        # Note that Bowtie reverse complements any sequence aligning to the reverse strand.  
        # and so read[3] is the 3'-end of minus strand reads 

        # Filter to get rid of reads of particular length. Or a particular strand.
        if (length < minlength or length > maxlength):
            continue

        mapped_reads += 1

        # 16 is the minus strand, 0 is the plus strand
        if (read[1] == '16'):
            start = startp
            density_minus[chrom][start] += 1 

        if (read[1] == '0'):
            start = startp + length - 1
            density_plus[chrom][start] += 1 
 
    density_plus[sequence.id]  = [float(i) * 1000000 / float(mapped_reads) for i in density_plus[sequence.id]]
    density_minus[sequence.id] = [float(i) * 1000000 / float(mapped_reads) for i in density_minus[sequence.id]]
    
    if source == 'TruSeq':
        makePickle(density_minus,path_den+"_plus")
        countstowig(density_minus,path_wig+"_plus")

        makePickle(density_plus,path_den+"_minus")
        countstowig(density_plus,path_wig+"_minus")    
    else:
        makePickle(density_plus,path_den+"_plus")
        countstowig(density_plus,path_wig+"_plus")

        makePickle(density_minus,path_den+"_minus")
        countstowig(density_minus,path_wig+"_minus")        

    print "\t" + fname + ' ' + str(mapped_reads)

    
def run_density(files, minlength, maxlength, threads, paths_in, paths_out, source): 
    
    print "\n-----DENSITY-----"
    print '\nFiles to condense: ' + ', '.join(files)
    print "\n\tStarted density at " + str(datetime.now())
    
    arguments = []
    
    for fname in files:

        # make paths for density files
          
        path_d = paths_out['path_density'] 
        path_w = paths_out['path_wig']

        path_den = path_d + fname
        path_wig = path_w + fname
        path_sam = paths_out['path_chr'] + fname + '_match.SAM'
        path_gff = paths_in['path_gff']
                       
        if not os.path.exists(path_sam):
            print "ERROR: " + fname + " has no alignment file, has been removed from analysis"
            continue
       
        argument = [fname, path_sam, minlength, maxlength, path_wig, path_den, path_gff, source]
        arguments.append(argument)
    
    multiprocess(density_3, arguments, threads)
            
    
    print "\tFinished density at " + str(datetime.now())

    print "\tCOMPLETED DENSITY"
    

    
    