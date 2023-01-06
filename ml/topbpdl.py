# Dataloader
from pathlib import Path
import pickle

import pandas as pd
import numpy as np
import pytorch_lightning as pl
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader, random_split
import ph_status_check
import matplotlib.pyplot as plt
import seaborn as sns

class TopBP_DL_Dataset(Dataset):
    def __init__(self,
        ph_status_df: pd.DataFrame):

        # filters dataframe for only correctly finished entries
        self.df = ph_status_df.loc[(ph_status_df['attempted'] == True) & (ph_status_df['finished'] == True) & (ph_status_df['error'] == False)]
        self.df.reset_index(inplace=True)

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        id = self.df.at[idx, 'name']
        item_base_folder = self.df.at[idx, 'folder']

        alpha_2dcnn_file = f'{item_base_folder}/{id}_protein_feature_complex_alpha_2DCNN.npy'
        alpha_2dcnn_data = torch.tensor(np.load(alpha_2dcnn_file), dtype=torch.float)

        distance_1dcnn_file = f'{item_base_folder}/{id}_protein_feature_complex_interaction_1DCNN.npy'
        distance_1dcnn_data = torch.tensor(np.load(distance_1dcnn_file), dtype=torch.float)
        distance_1dcnn_data = torch.swapaxes(distance_1dcnn_data, -1, -2)

        charge_1dcnn_file = f'{item_base_folder}/{id}_protein_feature_complex_electrostatics_1DCNN.npy'
        charge_1dcnn_data = torch.tensor(np.load(charge_1dcnn_file), dtype=torch.float)
        charge_1dcnn_data = torch.swapaxes(charge_1dcnn_data, -1, -2)

        affinity = self.df.at[idx, '-logKd/Ki']
        affinity = torch.tensor([affinity], dtype=torch.float)

        return alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data, affinity

class TopBP_DL_DataModule(pl.LightningDataModule):
    def __init__(self,
        ph_status_df : pd.DataFrame,
        batch_size = 8,
        num_workers = 12) -> None:

        super().__init__()
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.dataset = TopBP_DL_Dataset(ph_status_df)

    def setup(self, stage: str):
        # Assign train/val datasets for use in dataloaders
        train_size = int(0.6 * len(self.dataset))
        val_size = int(0.2 * len(self.dataset))
        test_size = int(0.1 * len(self.dataset))
        predict_size = len(self.dataset) - train_size - val_size - test_size

        self.train_set, self.val_set, self.test_set, self.predict_set = random_split(self.dataset, [train_size, val_size, test_size, predict_size])

    def train_dataloader(self):
        return DataLoader(self.train_set, batch_size=self.batch_size, num_workers=self.num_workers)

    def val_dataloader(self):
        return DataLoader(self.val_set, batch_size=self.batch_size, num_workers=self.num_workers)

    def test_dataloader(self):
        return DataLoader(self.test_set, batch_size=self.batch_size, num_workers=self.num_workers)

    def predict_dataloader(self):
        return DataLoader(self.predict_set, batch_size=self.batch_size, num_workers=self.num_workers)

class TopBP_DL_Model(nn.Module):
    def __init__(self) -> None:
        super(TopBP_DL_Model, self).__init__()

        self.alpha_2dcnn_head = nn.Sequential(
            nn.Conv2d(16, 64, 3), nn.ReLU(), nn.Conv2d(64, 64, 3), nn.ReLU(),
            nn.AvgPool2d(2), nn.Dropout(0.25),
            nn.Conv2d(64, 128, 3), nn.ReLU(), nn.Conv2d(128, 128, 3), nn.ReLU(),
            nn.AvgPool2d(2), nn.Dropout(0.25),
            nn.Flatten()
        )

        self.distance_1dcnn_head = nn.Sequential(
            nn.Conv1d(36, 128, 3), nn.ReLU(), nn.Conv1d(128, 128, 3), nn.ReLU(),
            nn.AvgPool1d(2), nn.Dropout(0.25),
            nn.Conv1d(128, 256, 3), nn.ReLU(), nn.Conv1d(256, 256, 3), nn.ReLU(),
            nn.AvgPool1d(2), nn.Dropout(0.25),
            nn.Flatten()
        )

        self.charge_1dcnn_head = nn.Sequential(
            nn.Conv1d(50, 128, 3), nn.ReLU(), nn.Conv1d(128, 128, 3), nn.ReLU(),
            nn.AvgPool1d(2), nn.Dropout(0.25),
            nn.Conv1d(128, 256, 3), nn.ReLU(), nn.Conv1d(256, 256, 3), nn.ReLU(),
            nn.AvgPool1d(2), nn.Dropout(0.25),
            nn.Flatten()
        )

        self.finale = nn.Sequential(
            nn.Linear(117888, 4096), nn.Tanh(), nn.Linear(4096, 4096), nn.Tanh(),
            nn.Linear(4096, 4096), nn.ReLU(), nn.Linear(4096, 4096), nn.ReLU(),
            nn.Dropout(0.5), nn.Linear(4096, 1)
        )



    def forward(self, alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data):
        # Input sizes:
        # torch.Size([16, 120, 128])
        # torch.Size([36, 200])
        # torch.Size([50, 100])

        alphahead_output = self.alpha_2dcnn_head(alpha_2dcnn_data)
        distancehead_output = self.distance_1dcnn_head(distance_1dcnn_data)
        chargehead_output = self.charge_1dcnn_head(charge_1dcnn_data)

        # print([o.shape for o in [alphahead_output, distancehead_output, chargehead_output]])
        y = torch.concat([alphahead_output, distancehead_output, chargehead_output], dim=1)
        # print(y.shape)
        y = self.finale(y)

        return y


class TNet_BP_Model(nn.Module):
    def __init__(self) -> None:
        super(TNet_BP_Model, self).__init__()

        self.alpha_2dcnn_head = nn.Sequential(
            nn.Conv2d(16, 64, 3), nn.ReLU(), nn.Conv2d(64, 64, 3), nn.ReLU(),
            nn.AvgPool2d(2), nn.Dropout(0.25),
            nn.Conv2d(64, 128, 3), nn.ReLU(), nn.Conv2d(128, 128, 3), nn.ReLU(),
            nn.AvgPool2d(2), nn.Dropout(0.25),
            nn.Flatten()
        )

        self.distance_1dcnn_head = nn.Sequential(
            nn.Conv1d(36, 128, 3), nn.ReLU(), nn.Conv1d(128, 128, 3), nn.ReLU(),
            nn.AvgPool1d(2), nn.Dropout(0.25),
            nn.Conv1d(128, 256, 3), nn.ReLU(), nn.Conv1d(256, 256, 3), nn.ReLU(),
            nn.AvgPool1d(2), nn.Dropout(0.25),
            nn.Flatten()
        )

        self.charge_1dcnn_head = nn.Sequential(
            nn.Conv1d(50, 128, 3), nn.ReLU(), nn.Conv1d(128, 128, 3), nn.ReLU(),
            nn.AvgPool1d(2), nn.Dropout(0.25),
            nn.Conv1d(128, 256, 3), nn.ReLU(), nn.Conv1d(256, 256, 3), nn.ReLU(),
            nn.AvgPool1d(2), nn.Dropout(0.25),
            nn.Flatten()
        )

        self.finale = nn.Sequential(
            nn.Linear(117888, 4096), nn.Tanh(), nn.Linear(4096, 4096), nn.Tanh(),
            nn.Linear(4096, 4096), nn.ReLU(), nn.Linear(4096, 4096), nn.ReLU(),
            nn.Dropout(0.5), nn.Linear(4096, 1)
        )



    def forward(self, alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data):
        # Input sizes:
        # torch.Size([16, 120, 128])
        # torch.Size([36, 200])
        # torch.Size([50, 100])

        alphahead_output = self.alpha_2dcnn_head(alpha_2dcnn_data)
        distancehead_output = self.distance_1dcnn_head(distance_1dcnn_data)
        chargehead_output = self.charge_1dcnn_head(charge_1dcnn_data)

        # print([o.shape for o in [alphahead_output, distancehead_output, chargehead_output]])
        y = torch.concat([alphahead_output, distancehead_output, chargehead_output], dim=1)
        # print(y.shape)
        y = self.finale(y)

        return y


class TopBP_DL_Module(pl.LightningModule):
    def __init__(self):
        super().__init__()
        self.criterion = nn.MSELoss()
        self.model = TopBP_DL_Model()

    def forward(self, alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data):
        return self.model(alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data)

    def on_fit_start(self):
        self.tensorboard = self.logger.experiment
        self.tensorboard.add_text('Model', str(self))

    def training_step(self, batch, batch_idx):
        alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data, y = batch
        y_hat = self(alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data)
        loss = self.criterion(y_hat, y)
        self.log("train_loss", loss, on_step=True, on_epoch=True, prog_bar=True, logger=True)
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=0.002)

    def validation_step(self, batch, batch_idx):
        alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data, y = batch
        y_hat = self(alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data)
        loss = self.criterion(y_hat, y)
        self.log("val_loss", loss)

# %%
def train():
    m = TopBP_DL_Module()
    PH_STATUS_CSV_FILE = 'ph_status.csv'
    ph_status_df = ph_status_check.get_df(PH_STATUS_CSV_FILE)
    datamodule = TopBP_DL_DataModule(ph_status_df, batch_size=2)
    trainer = pl.Trainer(max_epochs=10, resume_from_checkpoint='lightning_logs/version_1/checkpoints/epoch=4-step=795.ckpt')
    trainer.fit(m, datamodule=datamodule)

def predict():
    PH_STATUS_CSV_FILE = 'ph_status.csv'
    ph_status_df = ph_status_check.get_df(PH_STATUS_CSV_FILE)
    ds = TopBP_DL_Dataset(ph_status_df)
    module = TopBP_DL_Module()
    module = module.load_from_checkpoint('lightning_logs/version_2/checkpoints/epoch=9-step=1590.ckpt')


    predicted = []
    actual = []
    for i in range(100):
    # for i in range(len(ds)):
        alpha_2dcnn_data, distance_1dcnn_data, charge_1dcnn_data, y = ds[i]
        y_hat = module(alpha_2dcnn_data[None, :, :], distance_1dcnn_data[None, :, :], charge_1dcnn_data[None, :, :])[0][0].detach().cpu().numpy()
        predicted.append(y_hat)

        actual.append(y[0].detach().cpu().numpy())
        print(y_hat, y[0].detach().cpu().numpy())


    predicted = np.array(predicted)
    actual = np.array(actual)

    save_base_folder = Path('plots')
    save_base_folder.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(np.stack((predicted, actual), axis=-1), columns=['Predicted -logKd/Ki', 'Actual -logKd/Ki'])
    df.to_csv(save_base_folder / 'outputs.csv')

    mse = np.sum((np.array(predicted) - np.array(actual))**2) / len(predicted)

    fig = plt.figure()
    sns.scatterplot(data=df, x='Actual -logKd/Ki', y='Predicted -logKd/Ki')
    ax = fig.gca()
    ax.set_title(f'MSE: {mse}')
    imfile = save_base_folder / 'predictions.jpg'
    plt.savefig(imfile)

if __name__ == '__main__':
    # train()
    predict()


