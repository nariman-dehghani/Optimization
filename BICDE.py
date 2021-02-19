"""
Author: Nariman Dehghani

The MyEnv class is the engine of the IBCDE algorithm developed
by Dehghani et al. (2021).


BICDE algorithm is developed via the integration of 
a binary differential evolutionary (BDE) algorithm with 
the improved (μ+λ)-constrained differential evolution (ICDE). 


You can modify it based on your optimization problem.

Details on this algorithm can be find in the following book chapter:
    "Evolutionary Optimization for Resilience-Based Planning 
    for Power Distribution Networks"
"""
import numpy as np

class MyEnv:
    
    def __init__(self):
        """
        Initialize according to your problem
        """


    def step(self, population):
        """
        ATTENTION (SHOULD BE UPDATED):
        This method accpets the population and returns
        the fitness anf constraints. It should be updated based
        on the objective function and cosntrainst of the problem
        """
        return fitness, Constraints
    
    def select_criterion(self, ETA, Constraints_temp):
        diff = np.max(Constraints_temp) - np.min(np.max(Constraints_temp, axis = 0))
        if diff < ETA:
            criteria = 1
        else:
            criteria = 2
        return criteria

    def constraint_violation(self, criteria, Constraint_eval):
        if criteria == 1:
            Con_Vio = np.sum(Constraint_eval, axis = 1)
        else:
            max_constaints = np.max(Constraint_eval, axis = 0)
            Constraint_eval = Constraint_eval / max_constaints[None,:]
            Con_Vio = np.sum(Constraint_eval, axis = 1)
            Con_Vio = Con_Vio / max_constaints.size
        return Con_Vio 

    def is_pareto_efficient(self, costs, return_mask = True):
        """
        Find the pareto-efficient points
        :param costs: An (n_points, n_costs) array
        :param return_mask: True to return a mask
        :return: An array of indices of pareto-efficient points.
            If return_mask is True, this will be an (n_points, ) boolean array
            Otherwise it will be a (n_efficient_points, ) integer array of indices.
        """
        is_efficient = np.arange(costs.shape[0])
        n_points = costs.shape[0]
        next_point_index = 0  # Next index in the is_efficient array to search for
        while next_point_index<len(costs):
            nondominated_point_mask = np.any(costs<costs[next_point_index], axis=1)
            nondominated_point_mask[next_point_index] = True
            is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
            costs = costs[nondominated_point_mask]
            next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
        if return_mask:
            is_efficient_mask = np.zeros(n_points, dtype = bool)
            is_efficient_mask[is_efficient] = True
            return is_efficient_mask
        else:
            return is_efficient
    
    def ArATM(self, Mu, Lambda, H, fitness_H, Con_Vio_H, A, criteria, N_TimeSteps, N_Components):
        [Index_Infeasible] = np.where(Con_Vio_H > 0)
        [Index_feasible] = np.where(Con_Vio_H <= 0)
        NP = Mu + Lambda
        if len(Index_Infeasible) == 0:
            type_pop = 1 # feasible
            idx_fintness_sort = np.argsort(fitness_H)
            Final_Storage = idx_fintness_sort[0:Mu]
            
            Population_nx = {}
            fitness_pop_nx = np.zeros(len(Final_Storage))
            Con_Vio_pop_nx = np.zeros(len(Final_Storage))
            index_pop_nx = -1             
            for key in Final_Storage:
                index_pop_nx += 1
                Population_nx[(index_pop_nx)] = H.get(key)
                fitness_pop_nx[index_pop_nx] = fitness_H[key]
                Con_Vio_pop_nx[index_pop_nx] = Con_Vio_H[key]
            
            A_nx = A.copy() 
            
        elif len(Index_Infeasible) == NP:
            type_pop = 2 # infeasible
            if len(A) != 0:
                fitness_A = np.zeros(len(A),)
                Con_Vio_A = np.zeros(len(A),)
                for ind in range(len(A)):
                    A_ind = A.get(ind)
                    R_Matrix_A = A_ind.reshape((N_TimeSteps, N_Components))
                    fitness_A[ind], Constraints = self.step(R_Matrix_A)
                    Con_Vio_A[ind] = self.constraint_violation(criteria, [Constraints])

                randsize = np.random.randint(len(A))
                rand_indices_A = np.random.randint(len(A), size=randsize)
                index_H = len(H) 
                rand_indices_A = rand_indices_A.astype(int)
                for key in rand_indices_A:
                    H[(index_H)] = A.get(key)
                    fitness_H = np.append(fitness_H,fitness_A[key])
                    Con_Vio_H = np.append(Con_Vio_H,Con_Vio_A[key])
                    index_H += 1
            
            costs = np.concatenate((fitness_H, Con_Vio_H))
            costs = costs.reshape(2,fitness_H.size)
            costs = costs.T
            COST = costs.copy()
            
            while_indicator = -1
            popsize_next = 0
            while popsize_next < Mu:
                while_indicator += 1
                is_efficient = self.is_pareto_efficient(COST, return_mask = False)
                rank_vio = np.argsort(Con_Vio_H[is_efficient])
                First_half = np.array(np.ceil(rank_vio.size/2))
                First_half = First_half.astype(int)
                store_idx_temp = is_efficient[rank_vio[0:First_half]]
                if while_indicator == 0:
                    Storage = store_idx_temp.copy()
                else:
                    store_indiv = Not_Stored[store_idx_temp]
                    Storage = np.append(Storage,store_indiv)
                Not_Stored = np.array([i for i in range(costs.shape[0]) if i not in Storage])
                COST = costs[Not_Stored,:]
                popsize_next = len(Storage)
                
            if  popsize_next > Mu:
                temp_diff = popsize_next - Mu
                Final_Storage = Storage[0:-temp_diff]
            else:
                Final_Storage = Storage.copy()

            Not_Stored_Final = np.array([i for i in range(costs.shape[0]) if i not in Final_Storage])

            Population_nx = {}
            fitness_pop_nx = np.zeros(len(Final_Storage))
            Con_Vio_pop_nx = np.zeros(len(Final_Storage))
            index_pop_nx = -1                      
            for key in Final_Storage:
                index_pop_nx += 1
                Population_nx[(index_pop_nx)] = H.get(key)
                fitness_pop_nx[index_pop_nx] = fitness_H[key]
                Con_Vio_pop_nx[index_pop_nx] = Con_Vio_H[key]

            A_nx = {} 
            fitness_A_nx = np.zeros(len(Not_Stored_Final))
            Con_Vio_A_nx = np.zeros(len(Not_Stored_Final))
            index_A_nx = -1
            for key in Not_Stored_Final:
                index_A_nx += 1
                A_nx[(index_A_nx)] = H.get(key)
        else:
            type_pop = 3 # semi-feasible
            Z1 = {}
            fitness_Z1 = np.zeros(len(Index_feasible))
            Con_Vio_Z1 = np.zeros(len(Index_feasible))
            index_z1 = -1
            for key in Index_feasible:
                index_z1 += 1
                Z1[(index_z1)] = H.get(key)
                fitness_Z1[index_z1] = fitness_H[key]
                Con_Vio_Z1[index_z1] = Con_Vio_H[key]
            
            fitness_best_z1 = np.min(fitness_Z1)
            best_idx_z1 = np.argmin(fitness_Z1)
            best_z1 = Z1[(best_idx_z1)]
            
            fitness_worst_z1 = np.max(fitness_Z1)
            worst_idx_z1 = np.argmax(fitness_Z1)
            worst_z1 = Z1[(worst_idx_z1)]
            
            Phi = len(Index_feasible)/NP
            Z2 = {}
            fitness_Z2 = np.zeros(len(Index_Infeasible))
            Con_Vio_Z2 = np.zeros(len(Index_Infeasible))
            index_z2 = -1
            for key in Index_Infeasible:
                index_z2 += 1
                Z2[(index_z2)] = H.get(key)
                coverted_fit = np.array([Phi*fitness_best_z1 + (1-Phi)*fitness_worst_z1])
                fitness_Z2[index_z2] = np.max([fitness_H[key],coverted_fit])
                Con_Vio_Z2[index_z2] = Con_Vio_H[key]
            
            fitness_Z1_nor = (fitness_Z1 - np.min([np.min(fitness_Z1),np.min(fitness_Z2)]))/(np.max([np.max(fitness_Z1),np.max(fitness_Z2)]) - np.min([np.min(fitness_Z1),np.min(fitness_Z2)]))
            fitness_Z2_nor = (fitness_Z2 - np.min([np.min(fitness_Z1),np.min(fitness_Z2)]))/(np.max([np.max(fitness_Z1),np.max(fitness_Z2)]) - np.min([np.min(fitness_Z1),np.min(fitness_Z2)]))
            Con_Vio_Z1_nor = Con_Vio_Z1.copy()
            if criteria == 1:
                Con_Vio_Z2_nor = (Con_Vio_Z2 - np.min(Con_Vio_Z2))/(np.max(Con_Vio_Z2) - np.min(Con_Vio_Z2))
            else:
                Con_Vio_Z2_nor = Con_Vio_Z2.copy()
            
            fitness_Z1_final = fitness_Z1_nor + Con_Vio_Z1_nor
            fitness_Z2_final = fitness_Z2_nor + Con_Vio_Z2_nor
            
            fintness_final = np.zeros(NP,)
            fintness_final[Index_feasible] = fitness_Z1_final
            fintness_final[Index_Infeasible] = fitness_Z2_final
            idx_fintness_final_sort = np.argsort(fintness_final)
            Final_Storage = idx_fintness_final_sort[0:Mu]
            
            Population_nx = {}
            fitness_pop_nx = np.zeros(len(Final_Storage))
            Con_Vio_pop_nx = np.zeros(len(Final_Storage))
            index_pop_nx = -1             
            for key in Final_Storage:
                index_pop_nx += 1
                Population_nx[(index_pop_nx)] = H.get(key)
                fitness_pop_nx[index_pop_nx] = fitness_H[key]
                Con_Vio_pop_nx[index_pop_nx] = Con_Vio_H[key]
            
            A_nx = A.copy() 
            
        return Population_nx, fitness_pop_nx, Con_Vio_pop_nx, A_nx

    def BICDE(self, ParentPop, ParentSize, OffspringSize, N_RVs, F_prob, C_prob, criteria, N_TimeSteps, N_Components):
        Qt = {}
        fitness_Q = np.zeros(OffspringSize,)
        Con_Vio_Q = np.zeros(OffspringSize,)
        items = list(range(ParentSize))
        for ind in range(ParentSize):
            # mutation
            [r1,r2,r3] = np.random.choice(items, 3, replace=False)
            mutant_points_temp1 = np.random.rand(N_RVs) < F_prob
            mutant_points_temp2 = np.array(ParentPop[(r2)]) != np.array(ParentPop[(r3)])
            mutant_points = np.logical_and(mutant_points_temp1, mutant_points_temp2)
            mutant = np.where(mutant_points, 1-np.array(ParentPop[(r1)]), np.array(ParentPop[(r1)]))

            # Binomial Crossover
            cross_points = np.random.rand(N_RVs) <= C_prob
            if not np.any(cross_points):
                cross_points[np.random.randint(0, N_RVs)] = True
            trial = np.where(cross_points, mutant, np.array(ParentPop[(ind)]))
            R_Matrix_trial = trial.reshape((N_TimeSteps, N_Components))
            fitness_Q[ind], Constraints = self.step(R_Matrix_trial)
            Con_Vio_Q[ind] = self.constraint_violation(criteria, [Constraints])
            Qt[(ind)] = trial
            
        return Qt, fitness_Q, Con_Vio_Q