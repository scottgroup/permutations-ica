import numpy as np
import pandas as pd


class Dataset():

    def __init__(self, path, data_slice):
        self.path = path
        self.data = None
        self.data_slice = data_slice

        # Loading dataset from path
        self.load_data()

        # Keep subselection of the data
        self.slice_dataset()

        # Managing NaN
        self.transform_NaN(drop=True)

        # Transform data to CPM
        self.data2CPM()

        # Simple scaling
        self.scaling_data()

        # Flip back to genes as rows
        self.data = self.data.T

    def load_data(self):
        """ Load dataset from tsv """
        self.data = pd.read_csv(
            self.path, sep='\t',
            index_col=[0, 1, 2, 3], header=[0, 1, 2, 3, 4, 5],
            na_values=[''],
            dtype={
                'HGNC_ID': 'object',
                'symbol': 'object',
                'ncbi_id': 'object',
                'ensembl_id': 'object',
            }
        )

        # Transpose the data, set the genes as columns
        self.data = self.data.transpose()

    def slice_dataset(self):
        """ For the variables in data_slice, only keep specified values """

        slice_list = list()
        for idx in self.data.index.names:
            if idx in self.data_slice:
                slice_list.append(self.data_slice[idx])
            else:
                slice_list.append(slice(None))

        self.data = self.data.loc[tuple(slice_list),]

    def transform_NaN(self, drop):
        """
        Either drops all gene with at least one NaN, or fills them with 0
        """
        if drop:
            self.data.dropna(axis=1, inplace=True)
        elif not drop:
            self.data.fillna(value=0, inplace=True)

    def data2CPM(self):
        """ Transform each in CPM """
        sum_row = self.data.sum(axis=1)
        self.data = self.data.div(sum_row, axis=0).mul(1e6, axis=0)

    def scaling_data(self):
        """ Scaling data """

        pseudo_count = 13
        # self.data[c] = np.log10(self.data[c] + pseudo_count)
        # self.data = np.log10(self.data + pseudo_count)
        self.data = 2*np.sqrt(self.data + 3.0/8.0)
