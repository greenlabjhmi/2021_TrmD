'''

    Functions for analyzing ribosome profiling data:
        -- making lists of genes with the number of associated reads
        -- calculating pauses / enriched ribosome density at various sites
        
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






import csv
from Bio import Seq
import struct
import math
import random
import numpy as np

def bootstraperror(list):
    cycles=100
    boot=[0,0]

    waveofmeans = [0 for x in range(cycles)]
    datapoints = len(list)
    tempwave = [0 for x in range(datapoints)]

    for i in range(cycles):
        for j in range(datapoints):
            #Get a randompoint
            tempwave[j] = list[random.randint(0,datapoints-1)]
        waveofmeans[i] = sum(tempwave)/float(datapoints)

    meanofmeans = sum(waveofmeans)/float(cycles)
    
    dev = [x - meanofmeans for x in waveofmeans]
    dev2 = [x*x for x in dev]
    sdofmeans = math.sqrt(sum(dev2)/(float(cycles)-1)) # for the standard deviation of the population

    return [meanofmeans,sdofmeans]


def writelisttoexcel(inlist,filestring):
    import csv
    writer = csv.writer(open(filestring+".csv", "wb"))
    for element in inlist:
        writer.writerow(element)

        
#Function to write out any binary list to disk.
def writelistbinint(list,f):    # f is a file handle, i.e. f=open("text.bin","wb")
    for element in list:
        f.write(struct.pack("f",element))
    f.close()        
        
        
def writedicttoexcel(genelists,filestring):
    writer = csv.writer(open(filestring+".csv", "wb"),delimiter=',')
    if type(genelists)!=list:
        genelists=[genelists]       ### converting any input dict to a length 1 list.

    for genelist in genelists:     
        if(genelist.has_key("headers")):   # This gets header as 1st line in csv.
            headerrecord=[]
            headerrecord.append("headers")
            for field in genelist["headers"]:
                headerrecord.append(field)
            writer.writerow(headerrecord)
        for gene in genelist.keys():
            generecord=[]
            generecord.append(gene)

            # New Feb 2013, check to see if we have a list or a single value.
            if type(genelist[gene])==list:
                for field in genelist[gene]:
                    generecord.append(field)
            else:
                generecord.append(genelist[gene])
            if gene == "headers":      # Skip since we did this above.
                continue
            writer.writerow(generecord)

            
def readindict(f):
    previousgene = ""
    counter = 1
    filegen = csv.reader(f,delimiter=',')
    output = {}
    for gene in filegen:
        if gene[0] == previousgene:
            modgenename = gene[0] + "_" + str(counter)
            counter += 1
        else:
            modgenename = gene[0]
            counter = 1
        output[modgenename] = []
        for column in gene[1:]:
            output[modgenename].append(column)
    previousgene = gene[0]
    f.close()
    return output


# This is a function that gives counts and sequences for a given feature number
# This function takes as input a chrom and a list of counts (list with plus and minus tracks)
# output is a list of 2. The first element is the counts and second the sequence.
# It does not check for overlap with nearby genes
# but does require the features to be CDS and not pseudogenes
def givegene(feat_num, chrom, counts, shift):

    if chrom.features[feat_num].sub_features == []:
        return [-1,-1]
    elif chrom.features[feat_num].sub_features[0].type != 'CDS':
        return [-1,-1]
    elif chrom.features[feat_num].qualifiers.has_key('pseudo'):
        return [-1,-1]

    if(chrom.features[feat_num].strand == 1):    # positive, plus strand
        strand = 0
    else:
        strand = 1        # minus strand

    # from the GFF this returns the ends of the gene
    # but always in the low to high number direction from the chrom
    # so 'start' is always lower than 'end' even on - strand, but backwards
    start = chrom.features[feat_num].location.start.position
    end = chrom.features[feat_num].location.end.position

    finalseq = Seq.Seq('')
    finalcounts = []

    # assign seq and shifted counts for + strand
    if strand == 0:
        finalseq += chrom[start:end].seq
        finalcounts += counts[strand][start - shift : end - shift]

    # assign seq and shifted counts for - strand
    elif strand == 1:
        finalseq += chrom[start:end].seq
        finalseq = finalseq.reverse_complement()
        finalcounts += counts[strand][start + shift : end + shift]
        finalcounts.reverse()

    # Check to see if there were no entries:
    if finalcounts == [] or str(finalseq) == '':
        return [-1,-1]

    return [finalcounts, finalseq]


def makegenelist(counts, chrom, shift, thresh, totalreads):
    
    genelist = {} ## The list we will output.
    genelist["headers"] = ["alias","feat_num","rpkm", "rpc", "reads"]

    missedthresh = 0
    illegalgenes = 0
    genesinlist = 0
    
    feat_num = 0
    
    for feature in chrom.features:

        # Import sequence and counts. This is using the new givegene which takes shift as an input parameter.
        gg = givegene(feat_num, chrom, counts, shift)
        genesequence = gg[1]
        genecounts = gg[0]

        # Get rid of dubious genes, nongenes, genes with overlap of others.
        if (genesequence == -1 or genecounts == -1):
            feat_num += 1
            illegalgenes += 1
            continue

        # Compute rpkm for each.
        rpkm = float(1000)*sum(genecounts)/len(genecounts)

        # compute reads per codon
        rpm_per_codon = 3 * sum(genecounts) / len(genecounts)
        readfactor = totalreads / float(1000000)
        readspercodon = readfactor * rpm_per_codon
        if readspercodon < thresh:
            feat_num += 1
            missedthresh += 1
            continue

        # to get counts per gene, not rpm or rpkm for DESeq2
        countspergene = sum(genecounts) * readfactor
        countspergene = int(round(countspergene))

        if "Alias" in feature.qualifiers:
            alias = feature.qualifiers["Alias"][0]
        elif "Name" in feature.qualifiers: # For coli 
            alias=feature.qualifiers["Name"][0]
        else:
            alias = "NA"

        genelist[feature.id]=[]
        genelist[feature.id].append(alias)
        genelist[feature.id].append(feat_num)
        genelist[feature.id].append(rpkm)
        genelist[feature.id].append(readspercodon)
        genelist[feature.id].append(countspergene)

        feat_num += 1
        genesinlist += 1

    print "Genes below threshold = " + str(missedthresh)
    print "Genes dropped by givegene (overlap, undesirable features, etc.) = " + str(illegalgenes)
    print "Genes in list = " + str(genesinlist)

    return genelist


def motifavg_ab_wf(chrom, counts, outfilestring, motifsize, inframe, thresh, totalreads, shift, avgwin, scorewin):
  
    avglist=[]
    outlist=[]
    outlist.append(["motif", "total_hits", "count_avg", "oldpause", "error_OP", "count_OP"])
    
    #   motifdata, for each key, new numbering 2020
    #0 = avg counts at site, to be normalized by the count of instances included
    #1 = total count of every appearance of the motif
    #2 = count of motif instances included in raw average value
    #3 = old pause score 
    #4 = bootstrap error of pause score distribution
    #5 = count of instances in old pause score
    
    motifdata = motifavg_ab(outfilestring, chrom, counts, motifsize, inframe, \
                thresh, totalreads, shift, avgwin, scorewin)

    for mm in motifdata.keys() :
    
        # normalize oldpause score and compute error
        if motifdata[mm][5] > 0:
            motifdata[mm][4] = bootstraperror(motifdata[mm][3])[1]
            motifdata[mm][3] = sum(motifdata[mm][3]) / float(motifdata[mm][5])
        
        # normalize raw counts
        if motifdata[mm][2] > 0:
            for i in range(sum(avgwin)):
                motifdata[mm][0][i] /= float(motifdata[mm][2])
                        
        outlist.append([mm] + motifdata[mm][1:6])
        avglist += (motifdata[mm][0])
        
    # Write out list.
    writelisttoexcel(outlist,outfilestring)    
        
    # Write out avg file.
    favg = open(outfilestring + ".bin","wb")
    writelistbinint(avglist,favg)
    favg.close()
    

def motifavg_ab(outfilestring, chrom, counts, motifsize, inframe, \
                thresh, totalreads, shift, avgwin, scorewin):

    motifdata={}    

    if inframe == 2:
        motiflen = (int(motifsize)) * 3
    elif inframe <= 1:
        motiflen = int(motifsize)
    
    # oldpause scores do not include density that is too close to the ends
    # endcut[0] from start or [1] from end
    endcut = [20,20]    
    
    feat_num = 0
    included = 0
    
    for feature in chrom.features:
        gg = givegene(feat_num, chrom, counts, shift) # returns 0 = counts, 1 = seq

        if gg[1] == -1:      # For genes that are not allowed, pseudo or not CDS
            feat_num += 1
            continue

        # Use the mean with  3'-end density, too many genes lost with median = 0
        genemean = np.mean(gg[0])    

        # compute the reads_per_codon for the gene, to be compared with thresh
        genelen = len(gg[0])
        rpm_per_codon = 3 * sum(gg[0]) / genelen 
        readfactor = totalreads / float(1000000)
        reads_per_codon = readfactor * rpm_per_codon
        if reads_per_codon >= thresh:
            included += 1    

        # Now go through the gene looking for motifs.
        position = 0
        while position < (genelen - motiflen):    
            currentsequence = gg[1][position:position + motiflen]
            if inframe == 2:
                currentsequence = str(currentsequence.translate())
            else:
                currentsequence = str(currentsequence)
            
            # Check if entry exists for this motif, if not, make it
            if not motifdata.has_key(currentsequence):
                motifdata[currentsequence] = [[0 for x in range(sum(avgwin))], 0, 0, [], 0, 0]
            # increment total hits count by one
            motifdata[currentsequence][1] += 1
            

            if reads_per_codon >= thresh:    

                # 1. Add counts to binary output of equally weighted, averaged hits.
                # These are divided by the mean so that binary output matches the oldpause score
                if (position - avgwin[0]) >= endcut[0] and \
                            (position + avgwin[1] + endcut[1]) <= genelen:    
                    if genemean > 0:
                        for i in range(sum(avgwin)):
                            motifdata[currentsequence][0][i] += (
                                (gg[0][position - avgwin[0] + i]) / genemean)    
                            # output will have the first position of motif (w/ shift)
                            # at the first point after midpoint, 
                            # ie position 35 if avgwin ranges 0-69.
                        motifdata[currentsequence][2] += 1    
                        # Increment count of averaged hits.
            
                # 2: calculate old pause score if non-zero mean
                # numerator for pause = density from scorewin [0] to [1]
                if genemean > 0 and (position + scorewin[0]) >= endcut[0] and \
                                      (position + scorewin[1] + endcut[1]) <= genelen: 
                    motifdata[currentsequence][3].append(
                        sum(gg[0][position + scorewin[0]: position + scorewin[1]])
                                  / (genemean * (scorewin[1] - scorewin[0])))
                    motifdata[currentsequence][5] += 1    
                    # Increment number of scored hits.

            if (inframe == 1) or (inframe == 2):
                position += 3
            elif inframe == 0:
                position += 1
        
        feat_num += 1
        
    comments =  "This is the output of the motifavg_ab_wf function.\n"
    comments += "This function finds instances of all motifs of a given length\n"
    comments += "in the genome and creates equal weighted average plots and scores.\n"
    comments += "Inframe = 0 looks in all frames, 1 only in frame, 2 for aa.\n"

    # Write output file of comments.
    fc = open(outfilestring + "_output.txt", "w")
    fc.write(comments)
    fc.write("\n")
    fc.write("motif size = " + str(motifsize) + "\n")
    fc.write("inframe = " + str(inframe) + "\n")
    fc.write("thresh = " + str(thresh) + " reads per codon\n")
    fc.write("genes included: " + str(included) + "\n")
    fc.write("total genes: " + str(feat_num) + "\n")
    fc.write("shift = "+ str(shift) + "\n")
    fc.write("average window = " + str(avgwin) + "\n")
    fc.write("score window = " + str(scorewin) + "\n")
    fc.write("endcut is: " + str(endcut) + "\n")
    fc.write("total reads = " + str(totalreads) + "\n")
    fc.close()

    return motifdata




