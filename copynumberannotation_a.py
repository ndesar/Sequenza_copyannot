"""

HOW TO RUN THIS SCRIPT? 

On the command line, to enable writing uncovered genes (default), type:

$ python copynumberannotation.py -seg 05_144HH1_BONE.copynumber_calls.txt -gtf gencode_protein_linc_GTFV17.txt -og annotatedcopynumber_gene.txt -oe annotatedcopynumber_exon.txt -flag T

To disable writing uncovered genes, type:

$ python copynumberannotation.py -seg 05_144HH1_BONE.copynumber_calls.txt -gtf gencode_protein_linc_GTFV17.txt -og annotatedcopynumber_gene.txt -oe annotatedcopynumber_exon.txt -flag F

"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-seg', '--segfile', type=str, help="seg filename", default='new refcode/k20_l100_13-101H1_LN_R_ILIAC.copynumber_calls.txt')#'05_144HH1_BONE.copynumber_calls.txt')
parser.add_argument('-gtf', '--gtffile', type=str, help="reference gtf filename", default='new refcode/gencode.v36.basic.annotation.gtf')#'gencode_protein_linc_GTFV17.txt')
parser.add_argument('-og', '--outputfile_gene', type=str, help="gene output filename", default='annotatedcopynumber_gene.txt')
#parser.add_argument('-oe', '--outputfile_exon', type=str, help="exon output filename", default='annotatedcopynumber_exon.txt')
parser.add_argument('-flag', '--print_flag', type=str , help="print uncovered genes?", default='T')

def checkrange(range_gtffile, range_segfile):
    ''' Function to check whether the range1 falls in the desired range2 '''
    
    gtf0, gtf1 = int(range_gtffile[0]), int(range_gtffile[1])
    seg0, seg1 = int(range_segfile[0]), int(range_segfile[1])
    
    # exact match
    if gtf0 >= seg0 and gtf1 <= seg1:
        return 1
    # partial match
    elif gtf0 < seg0 and gtf1 > seg0 and gtf1 <= seg1:
        return 0.3
    elif gtf0 >= seg0 and gtf0 <= seg1 and gtf1 > seg1:
        return 0.6
    # no match
    else:
        return 0
    
def main_func(fname_segfile,fname_gtffile,outputfilename_gene,print_uncovered_gene):
    # list of line numbers having 'gene'
    gene_indices = [] 
    
    with open(fname_gtffile) as f:    
        # read all lines
        lines_gtffile = f.readlines()   
        
        # loop over each line 
        for count, line in enumerate(lines_gtffile):
            
            if 'gene' in line.split('\t')[2]:
                # get the line numbers of the lines with 'gene'
                gene_indices.append(count)        
                
    f.close()
    
    with open(fname_segfile) as f:
        
        # read all lines
        lines_segfile = f.readlines()    
        
    f.close()
    
    outputfile_gene = open(outputfilename_gene, 'w')
    #outputfile_exon = open(outputfilename_exon, 'w')
        
    # loop over lines with 'gene' only
    for count, idx_gtffile in enumerate(gene_indices):
        
        print(count)
        # get the individual column entries of each line
        line_gtffile = lines_gtffile[idx_gtffile].rstrip().split('\t')
        
        # flag variable to check uncovered gene
        flag = True
            
        # loop over each line of the seg file
        for idx_segfile in range(1,len(lines_segfile)):
            
            # get the individual column entries of each line
            line_segfile = lines_segfile[idx_segfile].rstrip().split(' ')

            # if the chr number matches
            if line_gtffile[0] == line_segfile[1].strip("\""):
                
                # check the ranges of the two entries
                value = checkrange((line_gtffile[3],line_gtffile[4]),(line_segfile[2],line_segfile[3]))
                
                # for perfect match
                if value == 1:
                    
                    flag = False
                    
                    # print the lines in the output file

                    #print(line_gtffile)
                    genes = line_gtffile[-1].split(';')[:3]
                    genes[0] = genes[0].replace('gene_id "','').rstrip('"')
                    genes[1] = genes[1].replace('gene_type "','').rstrip('"')
                    genes[2] = genes[2].replace('gene_name "','').rstrip('"')
                    line_gtffile_to_write = line_gtffile[:5] + genes
                    outputfile_gene.write('\t'.join(line_gtffile_to_write))

                    line_segfile_to_write = line_segfile[-4:]
                    outputfile_gene.write('\t'+ '\t'.join(line_segfile_to_write) + '\tfcov\n')                    
                    
                    # leave the loop since no need to check further
                    break
                
                # for no match
                elif value == 0.3:

                    flag = False

                    genes = line_gtffile[-1].split(';')[:3]
                    genes[0] = genes[0].replace('gene_id "','').rstrip('"')
                    genes[1] = genes[1].replace('gene_type "','').rstrip('"')
                    genes[2] = genes[2].replace('gene_name "','').rstrip('"')
                    line_gtffile_to_write = line_gtffile[:5] + genes

                    line_gtffile_to_write[3] = line_segfile[2]

                    outputfile_gene.write('\t'.join(line_gtffile_to_write))

                    line_segfile_to_write = line_segfile[-4:]
                    outputfile_gene.write('\t'+ '\t'.join(line_segfile_to_write) + '\tpcov\n')
                    
                
                # for partial match
                elif value == 0.6:

                    flag = False

                    genes = line_gtffile[-1].split(';')[:3]
                    genes[0] = genes[0].replace('gene_id "','').rstrip('"')
                    genes[1] = genes[1].replace('gene_type "','').rstrip('"')
                    genes[2] = genes[2].replace('gene_name "','').rstrip('"')
                    line_gtffile_to_write = line_gtffile[:5] + genes

                    line_gtffile_to_write[4] = line_segfile[3]

                    outputfile_gene.write('\t'.join(line_gtffile_to_write))

                    line_segfile_to_write = line_segfile[-4:]
                    outputfile_gene.write('\t'+ '\t'.join(line_segfile_to_write) + '\tpcov\n')
                    
                    '''
                    partial_gene_text = 'partial_'
                   
                    # loop over the lines with 'exon'
                    for idx in range(idx_gtffile+1,gene_indices[count+1]):
                        
                        # get the individual column entries
                        line = lines_gtffile[idx].rstrip().split('\t')
                        
                        # check the range again
                        val = checkrange((line[3],line[4]),(line_segfile[2],line_segfile[3]))
                        
                        # for perfect match
                        if val == 1:
                            
                            # print the lines in the output file  
                            partial_gene_text += line[16].split(' ')[2] + ','
                            #outputfile_exon.write('covered\t'+lines_gtffile[idx][:-1])
                            #outputfile_exon.write(lines_segfile[idx_segfile].replace(' ','\t'))
                            
                            
                        elif val == 0.5:
                            
                            # print the lines in the output file
                            partial_gene_text += line[16].split(' ')[2] + ','
                            #outputfile_exon.write('partial\t'+lines_gtffile[idx][:-1])
                            #outputfile_exon.write(lines_segfile[idx_segfile].replace(' ','\t'))

                            if idx_segfile < len(lines_segfile)-1:
                                # Get the next line to check for partial match
                                GetNextLine = lines_segfile[idx_segfile+1].rstrip().split(' ')
                                # check range with the newly obtained line
                                val2 = checkrange((line[3],line[4]),(GetNextLine[2],GetNextLine[3]))
                                
                                if val2 == 0.5:
                                    # print the lines in the output file
                                    partial_gene_text += line[16].split(' ')[2] + ','
                                    #outputfile_exon.write('partial\t'+lines_gtffile[idx][:-1])
                                    #outputfile_exon.write(lines_segfile[idx_segfile+1].replace(' ','\t'))
                                           
                        else:
                            # print the lines in the output file
                            pass
                            #outputfile_exon.write('uncovered\t'+lines_gtffile[idx][:-1])
                            #outputfile_exon.write('\tuncovered'*(len(line_segfile))+'\n')
                        '''
                    
                    # print the lines in the output file  
                    # outputfile_gene.write(partial_gene_text.rstrip(','))
                    # outputfile_gene.write('\t'+lines_gtffile[idx_gtffile][:-1])
                    # outputfile_gene.write(lines_segfile[idx_segfile].replace(' ','\t'))

                    # go to outermost loop
                    #break

                # for no match
                else:
                    continue


        if flag == True and print_uncovered_gene == True:
            # print the lines in the output file 

            genes = line_gtffile[-1].split(';')[:3]
            genes[0] = genes[0].replace('gene_id "','').rstrip('"')
            genes[1] = genes[1].replace('gene_type "','').rstrip('"')
            genes[2] = genes[2].replace('gene_name "','').rstrip('"')
            line_gtffile_to_write = line_gtffile[:5] + genes
            outputfile_gene.write('\t'.join(line_gtffile_to_write)) 

            outputfile_gene.write('\t'+ '.\t'*4 + '\tucov\n')
                       
    outputfile_gene.close()
    #outputfile_exon.close()

if __name__ == '__main__':
    args = parser.parse_args()
    if args.print_flag == 'T':
        print_uncovered_gene = True
    else:
        print_uncovered_gene = False
    main_func(args.segfile,args.gtffile,args.outputfile_gene,print_uncovered_gene)    