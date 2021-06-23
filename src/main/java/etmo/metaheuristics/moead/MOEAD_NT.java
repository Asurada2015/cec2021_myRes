package etmo.metaheuristics.moead;

import Jama.Matrix;
import etmo.core.*;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import jmetal.metaheuristics.moead.Utils;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.Vector;

public class MOEAD_NT extends Algorithm {



    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public MOEAD_NT(ProblemSet problemSet) {
        super(problemSet);
        functionType_ = "_TCHE1";
    }

    private int populationSize_;
    /**
     * Stores the population
     */
    private SolutionSet population_;
    /**
     * Z vector (ideal point)
     */
    double[] z_;
    /**
     * Lambda vectors
     */
    // Vector<Vector<Double>> lambda_ ;
    double[][] lambda_;
    /**
     * T: neighbour size
     */
    int T_;
    /**
     * Neighborhood
     */
    int[][] neighborhood_;
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
    /**
     * Operators
     */
    Operator crossover_;
    Operator mutation_;

    String dataDirectory_;





    @Override
    public SolutionSet execute() throws JMException, ClassNotFoundException {
        int maxEvaluations;

        evaluations_ = 0;
        maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
        populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
        dataDirectory_ = this.getInputParameter("dataDirectory").toString();

//        多任务的个体
        population_ = new SolutionSet(populationSize_ * problemSet_.size());
//        indArray_ = new Solution[problemSet_.get(0).getNumberOfObjectives()];

        T_ = ((Integer) this.getInputParameter("T")).intValue();
        nr_ = ((Integer) this.getInputParameter("nr")).intValue();
        delta_ = ((Double) this.getInputParameter("delta")).doubleValue();

        neighborhood_ = new int[populationSize_][T_];
        //理想点不同
        z_ = new double[problemSet_.get(0).getNumberOfObjectives() * problemSet_.size()];

        // lambda_ = new Vector(problem_.getNumberOfObjectives()) ;
        lambda_ = new double[populationSize_][problemSet_.get(0).getNumberOfObjectives()];

        crossover_ = operators_.get("crossover"); // default: DE crossover
        mutation_ = operators_.get("mutation"); // default: polynomial mutation

        // STEP 1. Initialization
        // STEP 1.1. Compute euclidean distances between weight vectors and find
        // T
        initUniformWeight();
        // for (int i = 0; i < 300; i++)
        // System.out.println(lambda_[i][0] + " " + lambda_[i][1]) ;

        initNeighborhood();

        // STEP 1.2. Initialize population
        initPopulation();

        // STEP 1.3. Initialize z_
        initIdealPoint();



        // STEP 2. Update
        do {





            int tasknum = problemSet_.size();
            int[] permutation = new int[populationSize_ * tasknum];
            Utils.randomPermutation(permutation, populationSize_ * tasknum);

            for (int i = 0; i < populationSize_ * tasknum; i++) {
                int n = permutation[i]; // or int n = i;

                int task = n / populationSize_;
                int neiborN = n % populationSize_;

                int type;
                double rnd = etmo.util.PseudoRandom.randDouble();

                // STEP 2.1. Mating selection based on probability
                if (rnd < delta_) // if (rnd < realb)
                {
                    type = 1; // neighborhood
                } else {
                    type = 2; // whole population
                }
                Vector<Integer> p = new Vector<Integer>();
                matingSelection(p, neiborN, 2, type, task);

                // STEP 2.2. Reproduction
                Solution child;
                Solution[] parents = new Solution[3];

                parents[0] = population_.get(p.get(0));
                parents[1] = population_.get(p.get(1));

//                改成随机任务交叉效果更差
//                parents[0] = population_.get((p.get(0) + PseudoRandom.randInt(0, tasknum - 1) * populationSize_) % (tasknum * populationSize_));
//                parents[1] = population_.get((p.get(1) + PseudoRandom.randInt(0, tasknum - 1) * populationSize_) % (tasknum * populationSize_));

                parents[2] = population_.get(n);

                // Apply DE crossover
                child = (Solution) crossover_.execute(new Object[] { population_.get(n), parents });

                // Apply mutation
                mutation_.execute(child);

                // Evaluation
                problemSet_.get(task).evaluate(child);
                evaluations_++;

                // STEP 2.3. Repair. Not necessary

                // STEP 2.4. Update z_
                updateReference(child, task);

                // STEP 2.5. Update of solutions
                updateProblem(child, neiborN, type, task);


            } // for



            for (int i = 0; i < populationSize_ * tasknum; i++) {
                int n = permutation[i]; // or int n = i;

                int task = n / populationSize_;
                int neiborN = n % populationSize_;

                int type;
                double rnd = etmo.util.PseudoRandom.randDouble();

                // STEP 2.1. Mating selection based on probability
                if (rnd < delta_) // if (rnd < realb)
                {
                    type = 1; // neighborhood
                } else {
                    type = 2; // whole population
                }
                Vector<Integer> p = new Vector<Integer>();
                matingSelection(p, neiborN, 2, type, task);

                // STEP 2.2. Reproduction
                Solution child;
                Solution[] parents = new Solution[3];

                parents[0] = population_.get(p.get(0));
                parents[1] = population_.get(p.get(1));

//                改成随机任务交叉效果更差
//                parents[0] = population_.get((p.get(0) + PseudoRandom.randInt(0, tasknum - 1) * populationSize_) % (tasknum * populationSize_));
//                parents[1] = population_.get((p.get(1) + PseudoRandom.randInt(0, tasknum - 1) * populationSize_) % (tasknum * populationSize_));

                parents[2] = population_.get(n);

                // Apply DE crossover
                child = (Solution) crossover_.execute(new Object[] { population_.get(n), parents });

                // Apply mutation
                mutation_.execute(child);

//                加入迁移
                for (int t = 0; t < tasknum; t++){
                    if (t == task) continue;
//
                    double tranRand = etmo.util.PseudoRandom.randDouble();
//                    if (maxEvaluations % 1000 == 0){
                    if (tranRand <= 0.1){
                        problemSet_.get(t).evaluate(child);
                        evaluations_++;
                        updateReference(child, t);
                        updateProblem2(child, neiborN, type, t);
                    }
                }

            } // for
//
        } while (evaluations_ < maxEvaluations);

        return population_;
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
            size = neighborhood_[id].length;
        } else {
            size = populationSize_;
        }
        int[] perm = new int[size];

        Utils.randomPermutation(perm, size);

        for (int i = 0; i < size; i++) {
            int k;
            if (type == 1) {
                k = neighborhood_[id][perm[i]];
            } else {
                k = perm[i]; // calculate the values of objective function
                // regarding the current subproblem
            }
            double f1, f2;

//            邻域内随机选择权重向量,我觉得会有重复
            f1 = fitnessFunction(population_.get(k + taskId * populationSize_), lambda_[k], taskId);
            f2 = fitnessFunction(indiv, lambda_[k], taskId);

//            依次选择权重向量
//            f1 = fitnessFunction(population_.get(k), lambda_[id]);
//            f2 = fitnessFunction(indiv, lambda_[id]);

            if (f2 < f1) {
                population_.replace(k + taskId * populationSize_, new Solution(indiv));
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

        if (functionType_.equals("_TCHE1")) {
            double maxFun = -1.0e+30;

            int turn = problemSet_.get(0).getNumberOfObjectives() * taskId;
            for (int n = 0 + turn; n < problemSet_.get(0).getNumberOfObjectives() + turn; n++) {
                double diff = Math.abs(individual.getObjective(n) - z_[n]);

                double feval;
                if (lambda[n % problemSet_.get(0).getNumberOfObjectives()] == 0) {
                    feval = 0.0001 * diff;
                } else {
                    feval = diff * lambda[n % problemSet_.get(0).getNumberOfObjectives()];
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

    private void updateReference(Solution individual, int tasknum) {
        int turn = problemSet_.get(0).getNumberOfObjectives();
        for (int n = 0 + turn * tasknum; n < turn * (tasknum + 1); n++) {
            if (individual.getObjective(n) < z_[n]) {
                z_[n] = individual.getObjective(n);

//                indArray_[n] = individual;
            }
        }

    }

    private void updateProblem2(Solution indiv, int id, int type, int taskId) {
        // indiv: child solution
        // id: the id of current subproblem
        // type: update solutions in - neighborhood (1) or whole population
        // (otherwise)
        int size;
        int time;

        time = 0;


        size = populationSize_;
        int[] perm = new int[size];

        Utils.randomPermutation(perm, size);

        for (int i = 0; i < size; i++) {
            int k = perm[i];

            double f1, f2;

//            邻域内随机选择权重向量,我觉得会有重复
            f1 = fitnessFunction(population_.get(k + taskId * populationSize_), lambda_[k], taskId);
            f2 = fitnessFunction(indiv, lambda_[k], taskId);

//            依次选择权重向量
//            f1 = fitnessFunction(population_.get(k), lambda_[id]);
//            f2 = fitnessFunction(indiv, lambda_[id]);

            if (f2 < f1) {
                population_.replace(k + taskId * populationSize_, new Solution(indiv));
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

    private void matingSelection(Vector<Integer> list, int cid, int size, int type, int tasknum) {
        // list : the set of the indexes of selected mating parents
        // cid : the id of current subproblem
        // size : the number of selected mating parents
        // type : 1 - neighborhood; otherwise - whole population
        int ss;
        int r;
        int p;

        ss = neighborhood_[cid].length;
        while (list.size() < size) {
            if (type == 1) {
                r = PseudoRandom.randInt(0, ss - 1);
                p = neighborhood_[cid][r];
                // p = population[cid].table[r];
            } else {
                p = PseudoRandom.randInt(0 , populationSize_ - 1 );
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
                list.addElement(p + tasknum * populationSize_);
            }
        }

    }

    private void matingSelection2(Vector<Integer> list, int cid, int size, int type, int tasknum) {
        // list : the set of the indexes of selected mating parents
        // cid : the id of current subproblem
        // size : the number of selected mating parents
        // type : 1 - neighborhood; otherwise - whole population
        int ss;
        int r;
        int p;

        ss = neighborhood_[cid].length;
        while (list.size() < size) {
            if (type == 1) {
                r = PseudoRandom.randInt(0, ss - 1);
                p = neighborhood_[cid][r];
                // p = population[cid].table[r];
            } else {
                p = PseudoRandom.randInt(0 , populationSize_ - 1 );
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
                list.addElement(p + tasknum * populationSize_);
            }
        }

    }


    private void initIdealPoint() throws ClassNotFoundException, JMException {
        for (int i = 0; i < problemSet_.get(0).getNumberOfObjectives() * problemSet_.size(); i++) {
            z_[i] = 1.0e+30;
//            indArray_[i] = new Solution(problemSet_);
//            problemSet_.get(0).evaluate(indArray_[i]);
//            evaluations_++;
        } // for

        for (int i = 0; i < populationSize_ * problemSet_.size(); i++) {
            int tasknum = i / populationSize_;
            updateReference(population_.get(i), tasknum);
        } // for

    }

    private void initPopulation() throws ClassNotFoundException, JMException {
        for (int i = 0; i < populationSize_ * problemSet_.size(); i++) {
            Solution newSolution = new Solution(problemSet_);
            problemSet_.get(i / populationSize_).evaluate(newSolution);
            evaluations_++;
            population_.add(newSolution);
        } // for
    }

    private void initNeighborhood() {
        double[] x = new double[populationSize_];
        int[] idx = new int[populationSize_];

        for (int i = 0; i < populationSize_; i++) {
            // calculate the distances based on weight vectors
            for (int j = 0; j < populationSize_; j++) {
                x[j] = Utils.distVector(lambda_[i], lambda_[j]);
                // x[j] = dist_vector(population[i].namda,population[j].namda);
                idx[j] = j;
                // System.out.println("x["+j+"]: "+x[j]+ ". idx["+j+"]:
                // "+idx[j]) ;
            } // for

            // find 'niche' nearest neighboring subproblems
            Utils.minFastSort(x, idx, populationSize_, T_);
            // minfastsort(x,idx,population.size(),niche);

            System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
        } // for
    }

    private void initUniformWeight() {
        if ((problemSet_.get(0).getNumberOfObjectives() == 2) && (populationSize_ <= 300)) {
            for (int n = 0; n < populationSize_; n++) {
                double a = 1.0 * n / (populationSize_ - 1);
                lambda_[n][0] = a;
                lambda_[n][1] = 1 - a;
            } // for
        } // if
        else {
            String dataFileName;
            dataFileName = "W" + problemSet_.get(0).getNumberOfObjectives() + "D_" + populationSize_ + ".dat";

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
                        lambda_[i][j] = value;
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
}
