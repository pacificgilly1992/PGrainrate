import csv

hello = [['Me','You'],['293', '219'],['13','15']]
length = len(hello[0])

with open('test1.csv', 'wb') as testfile:
    csv_writer = csv.writer(testfile)
    for y in range(length):
        csv_writer.writerow([x[y] for x in hello])