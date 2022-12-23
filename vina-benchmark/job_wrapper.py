import argparse
import pandas as pd
import dispatch_jobs

def job(df, idx):
    from vinascpdb import run_on_folder
    folder = df.at[idx, 'folder']
    run_on_folder(folder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv')
    parser.add_argument('--idx')

    args = parser.parse_args()

    CSV_FILE = args.csv
    idx = args.idx

    df = dispatch_jobs.get_df(CSV_FILE)

    try:
        job(df, idx)
    except:
        df.at[idx, 'error'] = True
        df.to_csv(CSV_FILE)

    df.at[idx, 'finished'] = True
    df.to_csv(CSV_FILE)
