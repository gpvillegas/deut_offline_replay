#  pseudo code for testing writing to csv file

import  csv
import os.path

header1 = ['run_num', 'T1_rates [kHz]', 'kin_type'] 
header2 = ['h_eff', 'cpuLT', 'Comments']
header =  header1+header2

print(header)

# rows of information, except the user comment
row = "%i %.3f %s" % (3288, 6.87, 'heep_coin')
row2 = "%.3f %i " % (5.212, 6643)

# read user comment
#comment = "this is !&@# is }[\| my & crazy comment? "
comment = input("Please enter any relevant comments for this run: ")


# clean user comment out of weird characters or spaces and replace them with _
specialChars = "!@#$%^&*()+={[]}|\:;,<>?/\" "

for specialChar in specialChars:
        comment = comment.replace(specialChar, '_')

# convert row to list
my_list = row.split()
print('my_list:', my_list)

# other possible sub-lists 
my_list2 = row2.split()
print('my_list2:', my_list2)

total_list = my_list + my_list2
print('total_list:', total_list)

print('comment:', comment)
# append comment to total_list
total_list.append(comment)

print(total_list)


fname_path='UTILS_CAFE/runlist/cafe-2022_runlist.csv'

# check if run list exists, else create it and add a header
if os.path.isfile(fname_path):
    print (fname_path," exists !")

    with open(fname_path, "a") as f:
        wr = csv.writer(f,delimiter=",")
        wr.writerow(total_list)
        
else:
    print (fname_path," does NOT exist ! \n Will create it and add a header")
    
    with open(fname_path, "a") as f:
        wr = csv.writer(f,delimiter=",")
        wr.writerow(header)
        wr.writerow(total_list)
        

f.close()

