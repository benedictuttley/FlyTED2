# Create dataframe and upload to mysql - Example

import pandas as pd
import numpy as np
# Script to write one Gene file to the php database:
my_frame = pd.read_excel('1-5-07/1-5-07.xls');

# Drop unnecessary columns:
my_frame.drop(columns = ['user'])

print(my_frame['link to file'])
my_frame['link to file'] = my_frame['link to file'].str.replace('\\','/')
print(my_frame['link to file'])
my_frame['link to file'] = my_frame['link to file'].str.replace('../../../../Pictures from microscope/Liz','')
print(my_frame['link to file'])

# Save in csv format:
#my_frame.to_csv('/home/benedict/PycharmProjects/test/formatted.txt', sep=' ', index=False, header=True)



# Correct image file paths:
# Info needed:
# 1) Date
# 2) Image Name
# 3) Probe

# Trial Attempt:

probe = 'CG8564'
image_name = 'wt01.bmp'
date = '1.5.07'

import datetime

def string_together(probe, image_name, date):
    image_name = image_name.replace('.bmp', '')
    new_path = (probe + "_" + image_name + "_" + date + ".bmp")
    return new_path


def ok(row):    # Whole row has been passed through:
    probe = row['probe']    # Fetch the probe
    print("Probe" + probe)

    date = row['date']  # Fetch and format the date
    date = date.date()
    date = str(date)
    if not np.isnat(np.datetime64(date)):
        date = datetime.datetime.strptime(date, '%Y-%m-%d').strftime('%-d.%-m.%-y')
        print(date)

        old_path = row['link to file']
        print(old_path)

        counter = 0
        image_name = ""
        for ch in old_path:
            if counter == 3:
                image_name = image_name + ch
            if ch == '/':
                counter = counter + 1

        print(image_name)
        t = string_together(probe, image_name, date)
        print(t)

        root = old_path.replace(image_name,'')
        news = root + t
        print(news)

        print("")
        return row['link to file'].replace(row['link to file'],news)

my_frame['link to file'] = my_frame.apply(ok, axis=1)
print("THE TEST")
print(my_frame['link to file'])


#print(string_together(probe, image_name, date))

# Insert dataframe into existing mysql table:
import sqlalchemy

# Connect to mysql and inject dataframe:
database_username = 'pmauser'
database_password = '[*3G4vlSp7'
database_ip = '127.0.0.1'
database_name = 'FlyTed2'
database_connection = sqlalchemy.create_engine('mysql+mysqlconnector://{0}:{1}@{2}/{3}'.format(database_username, database_password, database_ip, database_name))
my_frame.to_sql(con=database_connection, name='Demo', if_exists='replace')

