import re
import argparse
import os
import datetime

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
    # read the file
    for i in range(1):
        fin = open('job', "r")
        str2 = fin.read()
        fin.close()
        # the dictionary has target_word:replacement_word pairs
        word_dic = {
        '$number': str(i),
        }
        # call the function and get the changed text
        str3 = replace_words(str2, word_dic)
 
        # write changed text back out
        strjob='array.%03d'%i 
        fout = open(strjob, "w")
        fout.write(str3)
        fout.close()
        os.system("qsub "+strjob)
