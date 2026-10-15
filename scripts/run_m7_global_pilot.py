from concurrent.futures import ProcessPoolExecutor
import run_m7_global_metabolism as m
if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=3) as pool:list(pool.map(m.run,m.read('configs/m7-global-pilot-strict.json')))
