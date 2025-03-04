import glob
import re
import sys


# txt file with bams
#PATH = 'bam_files.txt'
PATH = sys.argv[1]

# txt file with samples to exclude
#TO_EXCLUDE = 'bam_files_to_exclude.txt'
TO_EXCLUDE = sys.argv[2]
   

def readTxtColumn(filename, n_column):
    fh = open(filename, 'r')
    linie = fh.readlines()
    k = [ele.split()[n_column] for ele in linie]
    return k
    
def writeBamList(list_of_paths):
    wh = open('10_00_get_bams.txt', 'w')
    for ele in list_of_paths:
        wh.write(ele + '\n')
    wh.flush()
    wh.close()
    
def removeSamples(big_list, sample_names):
    new_list = []
    for ele in big_list:
        if (ele in sample_names) or (ele.split('/')[-1] in sample_names):
            continue
        else:
            new_list.append(ele)
    return new_list
    

    
if __name__ == '__main__':
    
    
    # Removing samples if given & saving to a new list
    all_files = readTxtColumn(PATH,0)

    if TO_EXCLUDE != 'None':
        samples_to_exclude = readTxtColumn(TO_EXCLUDE, 0)
        bams = removeSamples(all_files, samples_to_exclude)
        writeBamList(bams)
    else:
        bams = all_files
        writeBamList(bams)
