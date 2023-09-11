import numpy as np
from geneticalgorithm2 import geneticalgorithm2 as ga
import time
import os
import multiprocessing
import threading
import portalocker

def submit_batch_job(batch_file, job_index):
    print(f'job {job_index} is submitting')
    os.system(f'sbatch {batch_file}')


folder_locks = [False] * 3
class Test:

    def find_available_sim_path(self):
        while True:
            for i in range(3):
                folder_path = f'test{i}'
                if not folder_locks[i]:
                    return folder_path

            time.sleep(10)

    def f(self, X):

        import math
        a = X[0]
        b = X[1]
        c = X[2]
        s = 0
        for i in range(10000):
            s += math.sin(a*i) + math.sin(b*i) + math.cos(c*i)

        while True:
            for i in range(3):
                try:
                    test_folder = f'test{i}'
                    lock_file = f'{test_folder}/test.lock'
                    with open(lock_file,'w') as lock_f:
                        submit_batch_job('test.sh', i)
                        portalocker.lock(lock_f, portalocker.LOCK_EX)

                        with open(f'{test_folder}/test.txt','a+') as f:
                            f.write(f'{s}\n')
                        portalocker.unlock(lock_f)
                    return s

                except Exception as e:
                    # print(e)
                    continue



algorithm_param = {'max_num_iteration': 50,
                   'population_size':100,
                   'mutation_probability':0.1,
                   'elit_ratio': 0.01,
                   'parents_portion': 0.3,
                   'crossover_type':'uniform',
                   'mutation_type': 'uniform_by_center',
                   'selection_type': 'roulette',
                   'max_iteration_without_improv':None}   
    
varbound = np.array([[-10,10]]*3)
test = Test()
model = ga(function=test.f, dimension=3, 
    variable_type='real',           
    variable_boundaries=varbound, 
    algorithm_parameters = algorithm_param)

########

# start = time.time()
# model.run()
# end = time.time()
# print(end - start)

start = time.time()
model.run(set_function= ga.set_function_multiprocess(test.f, n_jobs = 6))
# Wall time: 31.7 s
end = time.time()
print(end - start)