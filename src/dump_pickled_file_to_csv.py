import pickle

folder_name = 'myco_human_example'

df = pickle.load('./' + folder_name + '/step2_df.p')
df.to_csv('./' + folder_name + '/step2_df.p')