import sys


def generate_sequence(reference, variations, x, y):
    seq = ''
    pos = x
    
    #print(variations)
    #print(x, end="\t")
    #print(y)
    i = x
    while(i <= y):
        #print(i)
        if i in variations:
            ref, alt = variations[i]
            #if reference[i-1:i+len(ref)-1] == ref:
            #print("Found Var")
            seq += alt
            #print("ALT")
            #print(seq)
            i += len(ref) 
            #pos += len(ref)
            #else:
            #    seq += reference[i]
            #    pos += 1
        else:
            seq += reference[i-1]
            #print("NoALT")
            #print(seq)
            i = i + 1
            #pos += 1
    #print(seq)
    return seq

def read_fasta(file_path):
    fasta_dict = {}
    header = ""
    seq = ""
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if header:
                    fasta_dict[header] = seq
                header = line.strip().split()[0][1:]
                seq = ""
            else:
                seq += line.strip()
    fasta_dict[header] = seq
    return fasta_dict


def compare_vcf_files(reference, file1, file2):
    def get_alternate_sequence(ref, alt):
        """Helper function to get the alternate sequence from a reference and an alternate"""
        alt_seq = alt
        if len(ref) > len(alt):
            alt_seq = alt + ref[len(alt):]
        elif len(ref) < len(alt):
            alt_seq = alt[len(ref):]
        return alt_seq
    
    
    vcf1 = []
    vcf2 = []
    with open(file1, 'r') as vcf1_file:
        for line in vcf1_file:
            if not line.startswith("#"):
                vcf1.append(line.strip().split("\t"))
    with open(file2, 'r') as vcf2_file:
        for line in vcf2_file:
            if not line.startswith("#"):
                vcf2.append(line.strip().split("\t"))
    for entry1 in vcf1:
        found = False
        for entry2 in vcf2:
            entry1_var = {}
            entry2_var = {}
            entry1_var[int(entry1[1])] = (entry1[3], entry1[4])
            entry2_var[int(entry2[1])] = (entry2[3], entry2[4])
            if int(entry1[1]) == 172901 and int(entry2[1]) == 172898:
                print(entry1_var)
                print(entry2_var)
                print(generate_sequence(reference[entry1[0]], entry1_var, int(entry1[1])-25, int(entry1[1])+25))
                print(generate_sequence(reference[entry1[0]], entry2_var, int(entry1[1])-25, int(entry1[1])+25))
            if entry1[:5] == entry2[:5] or generate_sequence(reference[entry1[0]], entry1_var, int(entry1[1])-25, int(entry1[1])+25) == generate_sequence(reference[entry1[0]], entry2_var, int(entry1[1])-25, int(entry1[1])+25):
                found = True
                break
        if not found:
            print("VCF entry not found in second file: ", entry1)
    for entry2 in vcf2:
        found = False
        for entry1 in vcf1:
            entry1_var = {}
            entry2_var = {}
            entry1_var[int(entry1[1])] = (entry1[3], entry1[4])
            entry2_var[int(entry2[1])] = (entry2[3], entry2[4])
            if entry2[:5] == entry1[:5] or generate_sequence(reference[entry2[0]], entry2_var, int(entry2[1])-25, int(entry2[1])+25) == generate_sequence(reference[entry1[0]], entry1_var, int(entry2[1])-25, int(entry2[1])+25):
                found = True
                break
        if not found:
            print("VCF entry not found in first file: ", entry2)
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py file1.vcf file2.vcf ref.fa")
        sys.exit(1)
    file1, file2, reference_filename = sys.argv[1:]
    reference = read_fasta(reference_filename)
    compare_vcf_files(reference, file1, file2)
    
