import re
import argparse
import os
import datetime

def writelog(inputtext,outputfile):
    if os.path.exists(outputfile):
        outfile = open(outputfile,'a')
    else:
        outfile = open(outputfile,'wb')
    outfile.write(inputtext)
    outfile.write("\n")
    outfile.write("This case was created "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    outfile.write("\n")
    outfile.write("\n")
    outfile.close()
                          
def replace_words(text, word_dic):
    """
    take a text and replace words that match a key in a dictionary with
    the associated value, return the changed text
    """
    rc = re.compile('|'.join(map(re.escape, word_dic)))
    def translate(match):
	return word_dic[match.group(0)]
    return rc.sub(translate, text)
 
     
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case")
    parser.add_argument("exper")
    parser.add_argument("nodes", type=int)
    parser.add_argument("hours", type=int)
    parser.add_argument ("-m", "--minutes", dest='minutes', default='0', type=int)
    args=parser.parse_args()
    test_file = "job_2017"
    # read the file
    fin = open(test_file, "r")
    str2 = fin.read()
    fin.close()
    cores=str(22*int(args.nodes))
    # the dictionary has target_word:replacement_word pairs
    word_dic = {
    '$case': args.case,
    '$exper': args.exper,
    '$hours': "%02d"%args.hours,
    '$minutes': "%02d"%args.minutes,
    '$nodes': str(args.nodes),
    '$cores': cores}
 
    # call the function and get the changed text
    str3 = replace_words(str2, word_dic)
 
    # write changed text back out
    fout = open(args.case+"/"+args.exper+"job", "w")
    fout.write(str3)
    fout.close()
    logtext=args.exper
    writelog(logtext,args.case+"/"+args.exper+"joblog")           
    writelog(logtext+'\n jobs='+args.case+"/"+args.exper+"job",'globallog')   
    os.system("cd "+args.case+";qsub -A n02-GENESIS "+args.exper+"job")
