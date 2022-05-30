import os.path
import csv
from csv import writer
from csv import reader

filename='test.csv'
file_exists = os.path.isfile(filename)

print(file_exists)

header = ['row1', 'row2', 'row3']
row = [3288, 'heep_coin', 3.63]

with open(filename, 'a+') as csvfile:
    header = ['row1', 'row2', 'row3']
    writer = csv.DictWriter(csvfile, delimiter=',', lineterminator='\n',fieldnames=header)

    if not file_exists:
        writer.writeheader()  # file doesn't exist yet, write a header

    writer.writerow()
        
#    writer = csv.writer(f)
#    writer.writerow(header)
