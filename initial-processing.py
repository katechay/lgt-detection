#VERSION 1: 4/16/18 9:02 PM COMPLETE
#Take a given fasta file, extract identifier lines
#Put data lines in a second, temporary, file to be read
#From that file, output a line of data RE frequencies
#VERSION 2: 4/21/18 5:03 PM COMPLETE
#extend from version 1, use and keep track of rotating window
#VERSION 3: 4/23/18 1:19 PM COMPLETE
#extend from version 2, compare with data outside of window
#VERSION 4:
#extend from version 3, identify outlier sequences

#Method to find instances of given subsequence in given data string
def get_subseq_diff(subseq, dnasequence, TotalDiff):
    count = 0
    output_file.write(subseq)
    subseq_len = len(subseq)
    for i in range(0,length):
        is_in = dnasequence[i:(i+subseq_len)] == subseq
        if is_in:
            count = count+1
    diff = Averages[subseq]-count
    TotalDiff = TotalDiff+abs(diff)
    output_file.write(" "+str(count)+"\t")
    return TotalDiff

#Method to get averages per window across genome
def get_ave_subseq(subseq, dnasequence):
    count = 0
    output_file.write(subseq)
    subseq_len = len(subseq)
    for i in range(0,length):
        is_in = dnasequence[i:(i+subseq_len)] == subseq
        if is_in:
            count = count+1
    count = count/len(dnasequence)
    count = int(1000 * count)
    Averages[subseq] = count
    output_file.write(" "+str(count)+"\t")

#Method to process one window of DNA
def process_window(sequence):
    TotalDiff = 0
    TotalDiff=get_subseq_diff("GC", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("GT", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("GA", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("CG", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("CT", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("CA", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("TG", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("TC", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("TA", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("AG", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("AT", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("AC", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("GG", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("TT", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("CC", sequence, TotalDiff)
    TotalDiff=get_subseq_diff("AA", sequence, TotalDiff)

    output_file.write("Var: "+str(TotalDiff)+"\n")

#Open Files
fasta_name =  input("Enter a valid input file: ")
fasta_file = open(fasta_name)
temp_output = open("temp.txt", "w")
output_file  = open("results.txt", "w")

#Convert into pure DNA data
print("Extracting Data...")
for line in fasta_file:
    #separate and end analysis of id lines
    is_id_line = line.startswith(">")
    if is_id_line:
        continue
    #Process Data Lines into temp_output
    temp_output.write(line)
temp_output.close()

#Open pure DNA data for reading and add summary of whole
temp_output_file = open("temp.txt")
totalsequence = temp_output_file.read()
length = len(totalsequence)
numwin=int(length/1000)+1
output_file.write("Total NT "+str(length)+ " Total Windows "+str(numwin)+"\n")
print("Summarizing Data...")
#Get averages across the whole sequence
output_file.write("Averages: \n")
Averages={}
get_ave_subseq("GC", totalsequence)
get_ave_subseq("GT", totalsequence)
get_ave_subseq("GA", totalsequence)
get_ave_subseq("CG", totalsequence)
get_ave_subseq("CT", totalsequence)
get_ave_subseq("CA", totalsequence)
get_ave_subseq("TG", totalsequence)
get_ave_subseq("TC", totalsequence)
get_ave_subseq("TA", totalsequence)
get_ave_subseq("AG", totalsequence)
get_ave_subseq("AT", totalsequence)
get_ave_subseq("AC", totalsequence)
get_ave_subseq("GG", totalsequence)
get_ave_subseq("TT", totalsequence)
get_ave_subseq("CC", totalsequence)
get_ave_subseq("AA", totalsequence)
output_file.write("\nWindows: \n")

#separate into windows for processing
#Average gene length cited at ~1000 bp, so that will be size of window
i = 0
windownum = 1
while((i+1000)<length):
    print("processing window: "+ str(windownum)+ "/" + str(numwin))
    window = totalsequence[i: i + 1000]
    process_window(window)
    windownum = windownum + 1
    i = i+1000
window = totalsequence[i:]
process_window(window)
temp_output_file.close()
output_file.close()
