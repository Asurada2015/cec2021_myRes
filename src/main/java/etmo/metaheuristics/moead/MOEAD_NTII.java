package etmo.metaheuristics.moead;

import etmo.core.*;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import jmetal.metaheuristics.moead.Utils;

import javax.sound.midi.Soundbank;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.*;

import static org.apache.commons.lang3.RandomUtils.nextDouble;

public class MOEAD_NTII extends Algorithm {

    private int populationSize_;
    /**
     * Stores the population
     */
    private SolutionSet[] population_;
    /**
     * Z vector (ideal point)
     */
    double[][] z_;
    /**
     * Lambda vectors
     */
    // Vector<Vector<Double>> lambda_ ;
    double[][][] lambda_;
    /**
     * T: neighbour size
     */
    int T_;
    /**
     * Neighborhood
     */
    int[][][] neighborhood_;
    /**
     * delta: probability that parent solutions are selected from neighbourhood
     */
    double delta_;
    /**
     * nr: maximal number of solutions replaced by each child solution
     */
    int nr_;
    //    Solution[] indArray_;
    String functionType_;
    int evaluations_;

    int taskNum;
    /**
     * Operators
     */

    int rate;

    int[] generation;

    int neiborType;

    ArrayList[] savedValues;

    double rmp;

    Operator crossover_;
    Operator mutation_;

    String dataDirectory_;




    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public MOEAD_NTII(ProblemSet problemSet) {
        super(problemSet);
        functionType_ = "_TCHE1";
    }

    @Override
    public SolutionSet execute() throws JMException, ClassNotFoundException {
        return null;
    }


    public SolutionSet[] execute2() throws JMException, ClassNotFoundException {
        int maxEvaluations;
        taskNum = problemSet_.size();

        rate = 30;

        evaluations_ = 0;
        maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
        populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
        dataDirectory_ = this.getInputParameter("dataDirectory").toString();
        T_ = ((Integer) this.getInputParameter("T")).intValue();
        nr_ = ((Integer) this.getInputParameter("nr")).intValue();
        delta_ = ((Double) this.getInputParameter("delta")).doubleValue();
        rmp = ((Double) this.getInputParameter("rmp")).doubleValue();


        crossover_ = operators_.get("crossover"); // default: DE crossover
        mutation_ = operators_.get("mutation"); // default: polynomial mutation


        //        多任务的个体
        population_ = new SolutionSet[taskNum];
        neighborhood_ = new int[taskNum][][];
        z_ = new double[taskNum][];
        lambda_ = new double[taskNum][][];
        generation = new int[taskNum];
        Collection collection;
        savedValues = new ArrayList[taskNum];


        for (int i = 0; i < taskNum; i++){
            population_[i] = new SolutionSet(populationSize_);
            neighborhood_[i] = new int[populationSize_][T_];
            z_[i] = new double[problemSet_.get(i).getNumberOfObjectives()];
            lambda_[i] = new double[populationSize_][problemSet_.get(i).getNumberOfObjectives()];
            generation[i] = 0;
            savedValues[i] = new ArrayList<>(populationSize_);

            initUniformWeight(i);
            // for (int i = 0; i < 300; i++)
            // System.out.println(lambda_[i][0] + " " + lambda_[i][1]) ;

            initNeighborhood(i);

            // STEP 1.2. Initialize population
            initPopulation(i);

            // STEP 1.3. Initialize z_
            initIdealPoint(i);
        }


        do {
            for (int taskId = 0; taskId < taskNum; taskId++){


                List<Integer> order = tourSelection(10, taskId);

                for (int i = 0; i < order.size(); i++){
                    int subProblemId = order.get(i);
                    int n = subProblemId; // or int n = i;
                    // int n = i ; // or int n = i;
                    double rnd = PseudoRandom.randDouble();

                    // STEP 2.1. Mating selection based on probability
                    if (rnd < delta_) // if (rnd < realb)
                    {
                        neiborType = 1; // neighborhood
                    } else {
                        neiborType = 2; // whole population
                    }
                    Vector<Integer> p = new Vector<Integer>();
                    matingSelection(p, n, 2, neiborType, taskId);

                    // STEP 2.2. Reproduction
                    Solution child;
                    Solution[] parents = new Solution[3];

                    parents[0] = population_[taskId].get(p.get(0));
                    parents[1] = population_[taskId].get(p.get(1));
                    parents[2] = population_[taskId].get(n);

                    // Apply DE crossover
                    child = (Solution) crossover_.execute(new Object[] { population_[taskId].get(n), parents });


                    // Apply mutation
                    mutation_.execute(child);

                    int doTask = chooseTack(taskId);

                    // Evaluation
                    problemSet_.get(doTask).evaluate(child);

//                    迁移不进行评价
//                    if (doTask == taskId){
//                        evaluations_++;
//                    }
                    evaluations_++;


                    // STEP 2.3. Repair. Not necessary

                    // STEP 2.4. Update z_
                    updateReference(child, doTask);

                    // STEP 2.5. Update of solutions
                    updateProblem(child, n, neiborType, doTask);
                    }

                    generation[taskId]++;
                if (rate != 0 && generation[taskId] % rate == 0){
                    updateUtility(taskId);
                }

            }


//            迁移
//            for (int taskId = 0; taskId < taskNum; taskId++){
//                    //迁移
//                    for (int t = 0; t < taskNum; t++){
//                        if (t == taskId) continue;
//
//                        int[] permutation = new int[populationSize_];
//                        Utils.randomPermutation(permutation, populationSize_);
//
//                        for (int i = 0; i < populationSize_; i++) {
//                            int n = permutation[i]; // or int n = i;
//                            // int n = i ; // or int n = i;
//                            int type = 2; // whole population
//
//                            Vector<Integer> p = new Vector<Integer>();
//                            matingSelection(p, n, 2, type, taskId);
//
//                            // STEP 2.2. Reproduction
//                            Solution child;
//                            Solution[] parents = new Solution[3];
//
//                            parents[0] = population_[taskId].get(p.get(0));
//                            parents[1] = population_[taskId].get(p.get(1));
//                            parents[2] = population_[taskId].get(n);
//
//                            // Apply DE crossover
//                            child = (Solution) crossover_.execute(new Object[] { population_[taskId].get(n), parents });
//
//                            // Apply mutation
//                            mutation_.execute(child);
//
////
//                        double tranRand = PseudoRandom.randDouble();
//                        if (tranRand <= 0.1){
//                            problemSet_.get(t).evaluate(child);
//                            evaluations_++;
//                            updateReference(child, t);
//                            updateProblem(child, n, type, t);
//                        }
//                    }
//                } // for
//
//            }//for



        } while ( evaluations_ < maxEvaluations);

        return population_;
    }

    private int chooseTack(int taskId) {
        double rd = PseudoRandom.randDouble();

        if (rd < rmp){
            neiborType = 2;
            int id = PseudoRandom.randInt(0, taskNum - 1);
            while (id == taskId){
                id = PseudoRandom.randInt(0, taskNum - 1);
            }
            return id;
        }
        return taskId;
    }

    private void updateUtility(int taskId) {
        double f1, f2, uti, delta;
        for (int n = 0; n < populationSize_; n++) {

            f1 = fitnessFunction(population_[taskId].get(n), lambda_[taskId][n], taskId);
            f2 = fitnessFunction((Solution)savedValues[taskId].get(n), lambda_[taskId][n], taskId);

            delta = (f2 - f1) / f2;
            if (delta > 0.001) {
                population_[taskId].get(n).setUtility(1.0);
            } else {
//                uti = (0.95 + (0.05 * delta / 0.001)) * utility[n];
                uti = (0.95 + (0.05 * delta / 0.001)) * ((double) population_[taskId].get(n).getUtility());
                population_[taskId].get(n).setUtility(uti < 1.0 ? uti : 1.0);
            }

            savedValues[taskId].set(n, new Solution(population_[taskId].get(n)));
        }


    }

    private void initUniformWeight(int taskId) {
        if ((problemSet_.get(taskId).getNumberOfObjectives() == 2) && (populationSize_ <= 300)) {
            for (int n = 0; n < populationSize_; n++) {
                double a = 1.0 * n / (populationSize_ - 1);
                lambda_[taskId][n][0] = a;
                lambda_[taskId][n][1] = 1 - a;
            } // for
        } // if
        else {
            String dataFileName;
            dataFileName = "W" + problemSet_.get(taskId).getNumberOfObjectives() + "D_" + populationSize_ + ".dat";

            try {
                // Open the file
                FileInputStream fis = new FileInputStream(dataDirectory_ + "/" + dataFileName);
                InputStreamReader isr = new InputStreamReader(fis);
                BufferedReader br = new BufferedReader(isr);

                int numberOfObjectives = 0;
                int i = 0;
                int j = 0;
                String aux = br.readLine();
                while (aux != null) {
                    StringTokenizer st = new StringTokenizer(aux);
                    j = 0;
                    numberOfObjectives = st.countTokens();
                    while (st.hasMoreTokens()) {
                        double value = (new Double(st.nextToken())).doubleValue();
                        lambda_[taskId][i][j] = value;
                        // System.out.println("lambda["+i+","+j+"] = " + value)
                        // ;
                        j++;
                    }
                    aux = br.readLine();
                    i++;
                }
                br.close();
            } catch (Exception e) {
                System.out.println(
                        "initUniformWeight: failed when reading for file: " + dataDirectory_ + "/" + dataFileName);
                e.printStackTrace();
            }
        } // else

    }

    private void initNeighborhood(int taskId) {
        double[] x = new double[populationSize_];
        int[] idx = new int[populationSize_];

        for (int i = 0; i < populationSize_; i++) {
            // calculate the distances based on weight vectors
            for (int j = 0; j < populationSize_; j++) {
                x[j] = Utils.distVector(lambda_[taskId][i], lambda_[taskId][j]);
                // x[j] = dist_vector(population[i].namda,population[j].namda);
                idx[j] = j;
                // System.out.println("x["+j+"]: "+x[j]+ ". idx["+j+"]:
                // "+idx[j]) ;
            } // for

            // find 'niche' nearest neighboring subproblems
            Utils.minFastSort(x, idx, populationSize_, T_);
            // minfastsort(x,idx,population.size(),niche);

            System.arraycopy(idx, 0, neighborhood_[taskId][i], 0, T_);
        } // for
    }

    private void initPopulation(int taskId) throws ClassNotFoundException, JMException {
        for (int i = 0; i < populationSize_; i++) {
            Solution newSolution = new Solution(problemSet_);
            newSolution.setUtility(1.0);
            problemSet_.get(taskId).evaluate(newSolution);
            evaluations_++;
            population_[taskId].add(newSolution);

            savedValues[taskId].add(new Solution(newSolution));

        } // for
    }

    private void initIdealPoint(int taskId) throws ClassNotFoundException, JMException {
        for (int i = 0; i < problemSet_.get(taskId).getNumberOfObjectives(); i++) {
            z_[taskId][i] = 1.0e+30;
        } // for
        for (int i = 0; i < populationSize_; i++) {
            updateReference(population_[taskId].get(i), taskId);
        } // for
    }

    private void updateReference(Solution individual, int tasknum) {
        int turn = problemSet_.get(0).getNumberOfObjectives();
        for (int n = 0 ; n < problemSet_.get(tasknum).getNumberOfObjectives(); n++) {
            if (individual.getObjective(n + turn * tasknum) < z_[tasknum][n]) {
                z_[tasknum][n] = individual.getObjective(n + turn * tasknum);
            }
        }

    }

    private void matingSelection(Vector<Integer> list, int cid, int size, int type, int taskId) {
        // list : the set of the indexes of selected mating parents
        // cid : the id of current subproblem
        // size : the number of selected mating parents
        // type : 1 - neighborhood; otherwise - whole population
        int ss;
        int r;
        int p;

        ss = neighborhood_[taskId][cid].length;
        while (list.size() < size) {
            if (type == 1) {
                r = PseudoRandom.randInt(0, ss - 1);
                p = neighborhood_[taskId][cid][r];
                // p = population[cid].table[r];
            } else {
                p = PseudoRandom.randInt(0, populationSize_ - 1);
            }
            boolean flag = true;
            for (int i = 0; i < list.size(); i++) {
                if (list.get(i) == p) // p is in the list
                {
                    flag = false;
                    break;
                }
            }

            // if (flag) list.push_back(p);
            if (flag) {
                list.addElement(p);
            }
        }

    }

    private void updateProblem(Solution indiv, int id, int type, int taskId) {
        // indiv: child solution
        // id: the id of current subproblem
        // type: update solutions in - neighborhood (1) or whole population
        // (otherwise)
        int size;
        int time;

        time = 0;

        if (type == 1) {
            size = neighborhood_[taskId][id].length;
        } else {
            size = populationSize_;
        }
        int[] perm = new int[size];

        Utils.randomPermutation(perm, size);

        for (int i = 0; i < size; i++) {
            int k;
            if (type == 1) {
                k = neighborhood_[taskId][id][perm[i]];
            } else {
                k = perm[i]; // calculate the values of objective function
                // regarding the current subproblem
            }
            double f1, f2;

//            邻域内随机选择权重向量,我觉得会有重复
            f1 = fitnessFunction(population_[taskId].get(k), lambda_[taskId][k],taskId);
            f2 = fitnessFunction(indiv, lambda_[taskId][k],taskId);


            if (f2 < f1) {
                population_[taskId].replace(k, new Solution(indiv));
                // population[k].indiv = indiv;
                time++;
            }
            // the maximal number of solutions updated is not allowed to exceed
            // 'limit'
            if (time >= nr_) {
                return;
            }
        }


    }

    private double fitnessFunction(Solution individual, double[] lambda, int taskId) {
        double fitness;
        fitness = 0.0;

        int turn = problemSet_.get(0).getNumberOfObjectives() * taskId;

        if (functionType_.equals("_TCHE1")) {
            double maxFun = -1.0e+30;

            for (int n = 0; n < problemSet_.get(taskId).getNumberOfObjectives(); n++) {
                double diff = Math.abs(individual.getObjective(n + turn) - z_[taskId][n]);

                double feval;
                if (lambda[n] == 0) {
                    feval = 0.0001 * diff;
                } else {
                    feval = diff * lambda[n];
                }
                if (feval > maxFun) {
                    maxFun = feval;
                }
            } // for

            fitness = maxFun;
        } // if
        else {
            System.out.println("MOEAD.fitnessFunction: unknown type " + functionType_);
            System.exit(-1);
        }
        return fitness;


    }

    private List<Integer> tourSelection(int depth, int taskId) {
        List<Integer> selected = new ArrayList<Integer>();
        List<Integer> candidate = new ArrayList<Integer>();

        for (int n = 0; n < populationSize_; n++) {
            candidate.add(n);
        }

        while (selected.size() < (int) (populationSize_ / 5.0)) {
            int best_idd = (int) (PseudoRandom.randDouble() * candidate.size());
            int i2;
            int best_sub = candidate.get(best_idd);
            int s2;
            for (int i = 1; i < depth; i++) {
                i2 = (int) (PseudoRandom.randDouble()* candidate.size());
                s2 = candidate.get(i2);
                if ((double) population_[taskId].get(s2).getUtility() >
                        (double) population_[taskId].get(best_idd).getUtility()) {
                    best_idd = i2;
                    best_sub = s2;
                }
            }
            selected.add(best_sub);
            candidate.remove(best_idd);
        }
        return selected;
    }

}
