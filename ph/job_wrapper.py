import argparse
import pandas as pd
import traceback
import dispatch_jobs
import redis

DB = dispatch_jobs.get_db()

def job(key):
    # The part where the job actually runs, given df and idx as input
    from homology import calculate
    d = DB.hgetall(key)
    calculate(d['folder'], d['name'] + '_ligand', d['name'] + '_protein')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--key')

    args = parser.parse_args()
    key = args.key

    print('key', key)

    try:
        job(key)
        print('job finished')
        d = DB.hgetall(key)
        d['finished'] = 'True'
        d['error'] = 'False'
        DB.hset(key, mapping=d)

    except Exception as err:
        print(Exception, err)
        print(traceback.format_exc())
        print('job error')

        d = DB.hgetall(key)
        d['finished'] = 'True'
        d['error'] = 'True'
        DB.hset(key, mapping=d)
