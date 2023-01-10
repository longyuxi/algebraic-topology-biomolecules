import argparse
import pandas as pd
import dispatch_jobs

def job(df, idx):
    # The part where the job actually runs, given df and idx as input
    from homology import calculate
    calculate(df.at[idx, 'folder'], df.at[idx, 'name'] + '_ligand', df.at[idx, 'name'] + '_protein')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv')
    parser.add_argument('--idx')

    args = parser.parse_args()

    CSV_FILE = args.csv
    idx = int(args.idx)

    print(f'csv: {CSV_FILE}')
    print(f'idx: {idx}')

    try:
        df = dispatch_jobs.get_df(CSV_FILE)
        job(df, idx)
    except:
        print('job error')
        df = dispatch_jobs.get_df(CSV_FILE)
        df.at[idx, 'error'] = True
        df.to_csv(CSV_FILE, index=False)

    df = dispatch_jobs.get_df(CSV_FILE)
    print('job finished')
    df.at[idx, 'finished'] = True
    df.to_csv(CSV_FILE, index=False)