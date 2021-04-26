package etmo.metaheuristics.moea_ac;

import Jama.Matrix;
import etmo.core.*;
import etmo.operators.crossover.EGG;
import etmo.operators.crossover.SBXCrossover;
import etmo.operators.selection.RandomSelection;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.comparators.ObjectiveComparator;
import etmo.util.ranking.NondominatedRanking;
import etmo.util.ranking.Ranking;
import jmetal.metaheuristics.moead.Utils;

import java.util.ArrayList;
import java.util.List;

//mating selection with the selection restricted in the same cluster + integrated learning in recombination step
public class MaOEA_ACT extends Algorithm{
    private SolutionSet[] population_;
    private SolutionSet[] offspringPopulation_;
    private SolutionSet[] union_;

    private int populationSize_;
    int generations_;
    int maxGenerations_;

    Operator selection_;
    Operator crossover_;
    Operator mutation_;
    Operator learning_;

    private double[][] zideal_; //ideal point
    private double[][] znadir_;//Nadir point
    double[][][] extremePoints_; // extreme points
    private double p = 1.0;

    public MaOEA_ACT(ProblemSet problem) {
        super(problem);
        //zideal_ = new double[problem.getNumberOfObjectives()];
        //znadir_ = new double[problem.getNumberOfObjectives()];
    }

    public SolutionSet[] execute2() throws JMException, ClassNotFoundException {

        /*
         * step1: Basic Setting of this Algorithm
         */
        baiscSetting();
        /*
         * step2: Initialize the Population
         */
        initPopulation();
        initIdealPoint();  // initialize the ideal point

        initNadirPoint();    // initialize the nadir point

        initExtremePoints(); // initialize the extreme points
        /*
         * Enter the main loop，into the process of evolution
         */
        while (generations_ < maxGenerations_) {
            /*
             * step3 and step4:Mating Selection and Recombination to generate Offspring Populations
             */
            generateOffspringPopulation();
            /*
             * step5:Environmental Selection
             */
            environmentalSelection(generations_);

            generations_++;
        }

//        Ranking ranking = new NondominatedRanking(population_);
//        return ranking.getSubfront(0);
        return  population_;
    }



    public SolutionSet execute() throws JMException, ClassNotFoundException {

        /*
         * step1: Basic Setting of this Algorithm
         */
        baiscSetting();
        /*
         * step2: Initialize the Population
         */
        initPopulation();
        initIdealPoint();  // initialize the ideal point

        initNadirPoint();    // initialize the nadir point

        initExtremePoints(); // initialize the extreme points
        /*
         * Enter the main loop，into the process of evolution
         */
        while (generations_ < maxGenerations_) {
            /*
             * step3 and step4:Mating Selection and Recombination to generate Offspring Populations
             */
            generateOffspringPopulation();
            /*
             * step5:Environmental Selection
             */
            environmentalSelection(generations_);

            generations_++;
        }

//        Ranking ranking = new NondominatedRanking(population_);
//        return ranking.getSubfront(0);
        return population_[0];
    }

    /*
     * step1: Basic Setting of this Algorithm
     */
    public void baiscSetting(){
        generations_ = 0;
        maxGenerations_ = ((Integer) this.getInputParameter("maxGenerations")).intValue();
        populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
        mutation_  = operators_.get("mutation");
        crossover_ = operators_.get("crossover");
        selection_ = operators_.get("selection");
        learning_ = operators_.get("learning");

    }

    /*
     * step2: Initialize the Population
     */
    public void initPopulation() throws JMException, ClassNotFoundException {
        population_ = new SolutionSet[problemSet_.size()];
        for (int i = 0; i < problemSet_.size(); i++){
            population_[i] = new SolutionSet(populationSize_);
        }

        for (int j = 0; j < problemSet_.size(); j++){
            for (int i = 0; i < populationSize_; i++) {
                Solution newSolution = new Solution(problemSet_);
                problemSet_.get(j).evaluate(newSolution);
                problemSet_.get(j).evaluateConstraints(newSolution);
                population_[j].add(newSolution);
            } // for
        }

    } // initPopulation

    /*
     * step3 and step4:Mating Selection and Recombination to generate Offspring Populations
     */
    public void generateOffspringPopulation() throws JMException{

        offspringPopulation_ = new SolutionSet[problemSet_.size()];

        for (int task = 0; task < problemSet_.size(); task++){
            offspringPopulation_[task] = new SolutionSet(populationSize_);
            if (crossover_.getClass() == SBXCrossover.class){
                Solution[] parents = new Solution[2];
                for (int i = 0; i < (populationSize_); i++) {
                    // obtain parents
                    parents = (Solution[]) selection_.execute(population_[task]);

                    Solution[] offSpring = (Solution[]) crossover_
                            .execute(parents);
                    mutation_.execute(offSpring[0]);
                    problemSet_.get(task).evaluate(offSpring[0]);
                    problemSet_.get(task).evaluateConstraints(offSpring[0]);
                    offspringPopulation_[task].add(offSpring[0]);
                } // for
            }
            else if (crossover_.getClass() == EGG.class){

                int[] permutation = new int[populationSize_];
                Utils.randomPermutation(permutation, populationSize_);
                etmo.util.Ranking ranking = new etmo.util.Ranking(population_[task]);
                SolutionSet front = ranking.getSubfront(0);

                front.Suppress();
                SolutionSet KP = null;
                if(front.size()> problemSet_.get(task).getNumberOfObjectives())
                    KP = findingKneePoint(front, task);
                if(KP == null){
                    KP = population_[task];
                }

                for (int i = 0; i < (populationSize_); i++) {
                    // obtain parents
                    int n = permutation[i];
                    // STEP 2.1. Mating selection
                    int r1,r2;
                    r1 = PseudoRandom.randInt(0, KP.size() - 1);
                    do {
                        r2 = PseudoRandom.randInt(0, KP.size() - 1);
                    } while (r2==r1);
                    // STEP 2.2. Reproduction
                    Solution child;
                    Solution[] parents = new Solution[3];

                    parents[1] = KP.get(r1);
                    parents[2] = KP.get(r2);
                    parents[0] = population_[task].get(n);
                    child = (Solution) crossover_.execute(parents);

                    mutation_.execute(child);

                    problemSet_.get(task).evaluate(child);
                    problemSet_.get(task).evaluateConstraints(child);
                    offspringPopulation_[task].add(child);
                } // for
            }

        }




    }



    /*
     * step5:Environmental Selection
     */
    public void environmentalSelection(int G){
        /*
         * step5.1:Combine the Population and the Offspring Population
         */
        union_ = new SolutionSet[problemSet_.size()];

        for (int task = 0; task < problemSet_.size(); task++){
            union_[task] = ((SolutionSet) population_[task]).union(offspringPopulation_[task]);
            /*
             * step5.2:Normalization the Combined Population
             */
            SolutionSet[] st = getStSolutionSet(union_[task], populationSize_);

//		St  = st[1]
//		ND  = st[0]

            //estimateIdealPoint(st[1]);
            //estimateNadirPoint(st[1]);
            updateIdealPoint(st[0], task);
            updateNadirPoint(st[0], task);
            normalizationObjective(st[1], task);
            if(G >= 0.0*maxGenerations_){
                //p = estimation_Curvature(st[1]);
                p = estimation_Curvature(st[0], task);
                //p = 0.5;
            }
            if(G % 50 == 0){
//			System.out.println("p = "+p);
            }
            computeDistanceToIdealPoint(st[1],p, task);
            population_[task].clear();


            List<SolutionSet> list = new <SolutionSet>ArrayList();
            for(int i=0;i<st[1].size();i++){
                SolutionSet sols = new SolutionSet();
                sols.add(st[1].get(i));
                list.add(sols);
            }
            if(list.size() < populationSize_){
                System.out.println("ListSize4 = "+list.size());
                System.exit(0);
            }

            list = new HierarchicalClustering(list).clusteringAnalysis(populationSize_, p);
            if(list.size() != populationSize_){
                System.out.println("ListSize1 = "+list.size());
                System.exit(0);
            }

            bestSolutionSelection(list, task);

        }


    }

    /*
     * Estimate the Ideal Point
     */
//    public void estimateIdealPoint(SolutionSet solutionSet){
//        for (int task = 0; task < problemSet_.size(); task++){
//            for(int i=0; i<problemSet_.get(task).getNumberOfObjectives();i++){
//                zideal_[task][i] = 1.0e+30;
//                for(int j=0; j<solutionSet.size();j++){
//                    if(solutionSet.get(j).getObjective(i) < zideal_[task][i]){
//                        zideal_[task][i] = solutionSet.get(j).getObjective(i);
//                    }//if
//                }//for
//            }//for
//
//        }
//
//    }

    /*
     * Estimate the Nadir Point
     */
//    public void estimateNadirPoint(SolutionSet solutionSet){
//        for(int i=0; i<problemSet_.get(task).getNumberOfObjectives();i++){
//            znadir_[i] = -1.0e+30;
//            for(int j=0; j<solutionSet.size();j++){
//                if(solutionSet.get(j).getObjective(i) > znadir_[i]){
//                    znadir_[i] = solutionSet.get(j).getObjective(i);
//                }//if
//            }//for
//        }//for
//    }

    void copyObjectiveValues(double[] array, Solution individual, int task) {
        int move = 0;
        for (int m = 0; m < task; m++){
            move += problemSet_.get(m).getNumberOfObjectives();
        }
        for (int i = 0; i < problemSet_.get(task).getNumberOfObjectives(); i++) {
            array[i] = individual.getObjective(i + move);
        }
    }
    double asfFunction(Solution sol, int j, int task) {

        int move = 0;
        for (int m = 0; m < task; m++){
            move += problemSet_.get(m).getNumberOfObjectives();
        }

        double max = Double.MIN_VALUE;
        double epsilon = 1.0E-6;

        int obj = problemSet_.get(task).getNumberOfObjectives();

        for (int i = 0; i < obj; i++) {

            double val = Math.abs((sol.getObjective(i + move) - zideal_[task][i])
                    / (znadir_[task][i] - zideal_[task][i]));

            if (j != i)
                val = val / epsilon;

            if (val > max)
                max = val;
        }

        return max;
    }

    double asfFunction(double[] ref, int j, int task) {
        double max = Double.MIN_VALUE;
        double epsilon = 1.0E-6;

        int obj = problemSet_.get(task).getNumberOfObjectives();

        for (int i = 0; i < obj; i++) {

            double val = Math.abs((ref[i] - zideal_[task][i])
                    / (znadir_[task][i] - zideal_[task][i]));


            if (j != i)
                val = val / epsilon;

            if (val > max)
                max = val;
        }

        return max;
    }

    void initIdealPoint() {
        zideal_ = new double[problemSet_.size()][];
        for (int i = 0; i < problemSet_.size(); i++){
            zideal_[i] = new double[problemSet_.get(i).getNumberOfObjectives()];
        }

        for (int num = 0; num < problemSet_.size(); num++){
            for (int j = 0; j < zideal_[num].length; j++) {
                int move = 0;
                for (int m = 0; m < num; m++){
                    move += zideal_[m].length;
                }
                zideal_[num][j] = Double.MAX_VALUE;
                for (int i = 0; i < populationSize_; i++) {
                    if (population_[num].get(i).getObjective(j + move) < zideal_[num][j])
                        zideal_[num][j] = population_[num].get(i).getObjective(j + move);
                }
            }
        }

    }

    void updateIdealPoint(SolutionSet pop, int task){

        int move = 0;
        for (int i = 0; i < task; i++){
            move += problemSet_.get(i).getNumberOfObjectives();
        }

        for (int j = 0; j < problemSet_.get(task).getNumberOfObjectives(); j++) {
            for (int i = 0; i < pop.size(); i++) {
                if (pop.get(i).getObjective(j + move) < zideal_[task][j])
                    zideal_[task][j] = pop.get(i).getObjective(j + move);
            }
        }
    }

    void initNadirPoint() {
        znadir_ = new double[problemSet_.size()][];
        for (int i = 0; i < problemSet_.size(); i++){
            znadir_[i] = new double[problemSet_.get(i).getNumberOfObjectives()];
        }

        for (int num = 0; num < problemSet_.size(); num++){
            int move = 0;
            for (int m = 0; m < num; m++){
                move += znadir_[m].length;
            }
            for (int j = 0; j < znadir_[num].length; j++) {
                znadir_[num][j] = Double.MIN_VALUE;
                for (int i = 0; i < populationSize_; i++) {
                    if (population_[num].get(i).getObjective(j + move) > znadir_[num][j])
                        znadir_[num][j] = population_[num].get(i).getObjective(j + move);
                }
            }
        }
    }

    public void initExtremePoints() {
        extremePoints_ = new double[problemSet_.size()][][];
        for (int i = 0; i < problemSet_.size(); i++){
            int obj = problemSet_.get(i).getNumberOfObjectives();
            extremePoints_[i] = new double[obj][obj];
            for (int k = 0; k < obj; k++){
                for (int j = 0; j < obj; j++){
                    extremePoints_[i][k][j] = 1.0e+30;
                }
            }

        }
    }

    void updateNadirPoint(SolutionSet pop, int task){

        updateExtremePoints(pop, task);


        int obj = problemSet_.get(task).getNumberOfObjectives();
        double[][] temp = new double[obj][obj];

        for (int i = 0; i < obj; i++) {
            for (int j = 0; j < obj; j++) {
                double val = extremePoints_[task][i][j] - zideal_[task][j];
                temp[i][j] = val;
            }
        }

        Matrix EX = new Matrix(temp);

        boolean sucess = true;

        if (EX.rank() == EX.getRowDimension()) {
            double[] u = new double[obj];
            for (int j = 0; j < obj; j++)
                u[j] = 1;

            Matrix UM = new Matrix(u, obj);

            Matrix AL = EX.inverse().times(UM);

            int j = 0;
            for (j = 0; j < obj; j++) {

                double aj = 1.0 / AL.get(j, 0) + zideal_[task][j];


                if ((aj > zideal_[task][j]) && (!Double.isInfinite(aj)) && (!Double.isNaN(aj)))
                    znadir_[task][j] = aj;
                else {
                    sucess = false;
                    break;
                }
            }
        }
        else
            sucess = false;


        if (!sucess){
            double zmax[] = computeMaxPoint(pop, task);
            for (int j = 0; j < obj; j++) {
                znadir_[task][j] = zmax[j];
            }
        }
    }

    public void updateExtremePoints(SolutionSet pop, int task){
        for (int i = 0; i < pop.size(); i++)
            updateExtremePoints(pop.get(i), task);
    }

    public void updateExtremePoints(Solution individual, int task){
        int obj = problemSet_.get(task).getNumberOfObjectives();
        for (int i = 0; i < obj; i++){
            double asf1 = asfFunction(individual, i, task);
            double asf2 = asfFunction(extremePoints_[task][i], i, task);

            if (asf1 < asf2){
                copyObjectiveValues(extremePoints_[task][i], individual, task);
            }
        }
    }

    double[] computeMaxPoint(SolutionSet pop, int task){

        int move = 0;
        for (int m = 0; m < task; m++){
            move += problemSet_.get(m).getNumberOfObjectives();
        }

        int obj = problemSet_.get(task).getNumberOfObjectives();
        double zmax[] = new double[obj];
        for (int j = 0; j < obj; j++) {
            zmax[j] = Double.MIN_VALUE;

            for (int i = 0; i < pop.size(); i++) {
                if (pop.get(i).getObjective(j + move) > zmax[j])
                    zmax[j] = pop.get(i).getObjective(j + move);
            }
        }
        return zmax;
    }

    /*
     * Normalization
     */
    public void normalizationObjective(SolutionSet solutionSet, int task){
        int move = 0;
        for (int m = 0; m < task; m++){
            move += problemSet_.get(m).getNumberOfObjectives();
        }
        for(int i=0; i<solutionSet.size(); i++){
            Solution sol = solutionSet.get(i);
            for(int j=0; j<problemSet_.get(task).getNumberOfObjectives(); j++){
                double val = 0.0;
                val = (sol.getObjective(j + move)-zideal_[task][j])/(znadir_[task][j]-zideal_[task][j]+1e-30);
                if (val == Double.NaN) System.out.println("set normalization is NaN!!");
                sol.setNormalizedObjective(j, val);
            }//for
        }//for
    }

    /*
     * Compute the Convergence Distance of each Solutions Which use the distance of
     * each solution to the Ideal Point
     */
    public void computeDistanceToIdealPoint(SolutionSet solutionSet, double p, int task){
        for(int i=0; i<solutionSet.size(); i++){
            Solution sol = solutionSet.get(i);
            double normDistance = 0.0;
            double sumValue = 0.0;
            double pValue = 0.0;
            for(int j=0; j<problemSet_.get(task).getNumberOfObjectives(); j++){
                normDistance += sol.getNormalizedObjective(j)*sol.getNormalizedObjective(j);
                sumValue +=  sol.getNormalizedObjective(j);
                pValue += Math.pow(sol.getNormalizedObjective(j), p);
            }
            normDistance = Math.sqrt(normDistance);
            pValue = Math.pow(pValue, 1.0/p);
            if (Double.isNaN((normDistance))) System.out.println("normDistance is NaN");
            sol.setDistanceToIdealPoint(normDistance);
//    		debug
            if (Double.isInfinite(sumValue)) System.out.println("sumValue = INFINITY");
            sol.setSumValue(sumValue);
            for(int j=0; j<problemSet_.get(task).getNumberOfObjectives(); j++){
                sol.setIthTranslatedObjective(j, sol.getNormalizedObjective(j)/(pValue + 1e-30));
            }
        }
    }

    public SolutionSet[] getStSolutionSet(SolutionSet ss,int size) {
        SolutionSet[] sets = new SolutionSet[2];
        Ranking ranking = new NondominatedRanking(ss);

        int remain = size;
        int index = 0;
        SolutionSet front = null;
        SolutionSet mgPopulation = new SolutionSet();
        front = ranking.getSubfront(index);
        sets[0] = front;
        while ((remain > 0) && (remain >= front.size())) {

            for (int k = 0; k < front.size(); k++) {
                mgPopulation.add(front.get(k));
            } // for

            // Decrement remain
            remain = remain - front.size();

            // Obtain the next front
            index++;
            if (remain > 0) {
                front = ranking.getSubfront(index);
            } // if
        }
        if (remain > 0) { // front contains individuals to insert
            for (int k = 0; k < front.size(); k++) {
                mgPopulation.add(front.get(k));
            }
        }

        sets[1] = mgPopulation;

        return sets;
    }

    public void bestSolutionSelection(List<SolutionSet> list, int task){
        boolean[] flag = new boolean[list.size()];
        SolutionSet solSet = new SolutionSet();
        for(int i=0;i<list.size();i++){
            for(int j=0;j<list.get(i).size();j++){
                list.get(i).get(j).setClusterID(i);
                solSet.add(list.get(i).get(j));
            }
        }
        SolutionSet cornerSet = new MaxSimilarBasedClustering(solSet, problemSet_.get(task).getNumberOfObjectives())
                .decomposePopulation(problemSet_.get(task).getNumberOfObjectives());
        for(int m=0;m<problemSet_.get(task).getNumberOfObjectives();m++){
            if(!flag[cornerSet.get(m).getClusterID()]){
                flag[cornerSet.get(m).getClusterID()] = true;
                population_[task].add(cornerSet.get(m));
            }/*else{
    		   System.out.println("重复添加了");
    		   System.exit(0);
    	   }*/
        }
        for(int k=0;k<list.size();k++){
            if(!flag[k]){
                SolutionSet sols = list.get(k);
                if(sols.size() == 0){
                    System.out.println("SolsSize2 = "+sols.size());
                    System.exit(0);
                }
                double minFitness = Double.MAX_VALUE;
                int minFitnessID = -1;
                if(sols.size() == 0){
                    System.out.println("size = 0!");
                    System.exit(0);
                }
                for(int j=0;j<sols.size();j++){
                    Solution sol2 = sols.get(j);
                    double fitness = sol2.getSumValue();

                    if (Double.isNaN(fitness)) System.out.println("fitness is NaN");
//   				System.out.printf("%d : %f\n",j, fitness);

                    if(minFitness > fitness){
                        minFitness = fitness;
                        minFitnessID = j;
                    }
                }//for
                population_[task].add(sols.get(minFitnessID));
            }
        }
    }

    public double estimation_Curvature(SolutionSet solutionSet, int task){
        SolutionSet solSet = solutionSet;
        double c = 1.0;
        int size = solSet.size();
        int numb = problemSet_.get(task).getNumberOfObjectives();

        double sum = 0.0;
        double mean = 0.0;
        double var = 0.0;
        double maxD = 0.0;

        SolutionSet sSet = new SolutionSet();
        for(int i=0;i<size;i++){
            Solution sol = solSet.get(i);
            sSet.add(sol);
        }

        for(int i=0;i<sSet.size();i++){
            Solution sol = sSet.get(i);
            for(int j=0; j<numb; j++){
                if(sol.getNormalizedObjective(j) > 1.0){
                    sSet.remove(i);
                    i--;
                    break;
                }
            }
        }
        solSet = sSet;
        size = sSet.size();

        double[] dis = new double[size];
        for(int i=0;i<size;i++){
            Solution sol = solSet.get(i);
            dis[i] = 0.0;
            for(int j=0; j<numb; j++){
                dis[i] += sol.getNormalizedObjective(j);
            }
            dis[i] = dis[i] - 1.0;
            dis[i] = dis[i]/Math.sqrt(numb);
            sum += dis[i];
        }

        mean = sum/size;

        sum = 0.0;
        for(int i=0;i<size;i++){
            sum += Math.pow(dis[i]-mean, 2);
        }

        if(size > 1){
            var = sum/(size-1);
        }else{
            var = sum/size;
        }
        var = Math.sqrt(var);
        double cv = var/Math.abs(mean);

        double dis1 = (Math.pow(numb, 1.0 - 1.0/0.5)-1.0)/(Math.sqrt(numb));
        double dis2 = (Math.pow(numb, 1.0 - 1.0/2.0)-1.0)/(Math.sqrt(numb));

        double alpha = 0.5 - problemSet_.get(task).getNumberOfObjectives()*0.02;

        if(mean < 0 && mean < alpha*dis1){
            c = 0.5;
        }else if(mean > 0 && mean > alpha*dis2){
            c = 2.0;
        }else{
            c = 1.0;
        }

        if(c != 1.0 && cv < 0.1){
            c = 1.0;
        }
        //c = 0.5;

        return c;
    }

    //   add from egg
    double r_;
    double t_;


    private SolutionSet findingKneePoint(SolutionSet pop, int task) {

        int move = 0;
        for (int m = 0; m < task; m++){
            move += problemSet_.get(m).getNumberOfObjectives();
        }

        Matrix Pop_M, ones_M, ExtremePoint_M, Hyperplane, Distance_M;
        SolutionSet KneePoint = new SolutionSet();
        SolutionSet ExtremePoint = new SolutionSet(problemSet_.get(task).getNumberOfObjectives());
        SolutionSet TempPop = new SolutionSet(pop.size());
        double [] f_max = new double [problemSet_.get(task).getNumberOfObjectives()];
        double [] f_min = new double [problemSet_.get(task).getNumberOfObjectives()];
        for (int i=0;i<problemSet_.get(task).getNumberOfObjectives();i++) {
            f_max[i]=pop.get(0).getObjective(i + move);
            f_min[i]=pop.get(0).getObjective(i + move);
        }
        for(int i=0;i<pop.size();i++) {

            Solution sol = pop.get(i);
            int start = problemSet_.get(task).getStartObjPos();
            int end = problemSet_.get(task).getEndObjPos();
            Solution newSolution = new Solution(end - start + 1);

            for (int l = start; l <= end; l++)
                newSolution.setObjective(l - start, sol.getObjective(l));
            TempPop.add(newSolution);

            for (int j=0;j<problemSet_.get(task).getNumberOfObjectives();j++) {
                if(pop.get(i).getObjective(j + move)>f_max[j]) f_max[j]=pop.get(i).getObjective(j + move);
                if(pop.get(i).getObjective(j + move)<f_min[j]) f_min[j]=pop.get(i).getObjective(j + move);
            }
        }
        ones_M = new Matrix(problemSet_.get(task).getNumberOfObjectives(),1,1.0);
        //1. find extrem solution
        for(int i=0;i<problemSet_.get(task).getNumberOfObjectives();i++){
            ObjectiveComparator oc = new ObjectiveComparator(i);
            int best = TempPop.indexBest(oc);
            ExtremePoint.add(TempPop.get(best));
            TempPop.remove(best);
        }
        for (int i=0;i<ExtremePoint.size();i++){
            KneePoint.add(ExtremePoint.get(i));
        }

        //2. calculate extreme hyperplane
        ExtremePoint_M = new Matrix(ExtremePoint.writeObjectivesToMatrix());
        if(!ExtremePoint_M.lu().isNonsingular()) return ExtremePoint;
        Hyperplane = ExtremePoint_M.inverse().times(ones_M);
        double normHyperplane = Hyperplane.norm2();
        //Calculate the distance between each solution in TempPop and L
        Pop_M = new Matrix(TempPop.writeObjectivesToMatrix());
        Distance_M = Pop_M.times(Hyperplane);

        //Calculate R point
        r_ = r_*Math.exp(-((1.0-t_/0.5)/(double)problemSet_.get(task).getNumberOfObjectives()));
        double[] R = new double[problemSet_.get(task).getNumberOfObjectives()];
        for(int i=0;i<problemSet_.get(task).getNumberOfObjectives();i++)
            R[i] = (f_max[i]-f_min[i])*r_;

        //Add KneePoint
        int remain = TempPop.size();
        boolean Choose[] = new boolean [Distance_M.getRowDimension()];
        for(int i=0;i<Distance_M.getRowDimension();i++) Choose[i]=false;
        while(remain>0){
            //find the largest distance in Distance_M that have contain in KneePoint
            double largest_dist = Double.NEGATIVE_INFINITY;
            int idx = -1;
            for(int i=0;i<Distance_M.getRowDimension();i++){
                if(!Choose[i]){// judge the ith individual selected/remove or not
                    double dist = -(Distance_M.get(i, 0)-1.0)/normHyperplane;
                    if(dist > largest_dist){
                        largest_dist = dist;
                        idx = i;
                    }
                }
            }
            if(idx!=-1){
                KneePoint.add(TempPop.get(idx));
                Solution currentKneeP = TempPop.get(idx);
                //remove neighbors
                for(int i=0;i<TempPop.size();i++){
                    if(!Choose[i]){
                        Solution s = TempPop.get(i);
                        boolean flag = true;
                        for(int j=0;j<problemSet_.get(task).getNumberOfObjectives();j++){
                            if(Math.abs(s.getObjective(j + move)-currentKneeP.getObjective(j + move))>R[j]) {flag = false;break;}
                        }
                        if(flag){
                            Choose[i]=true;remain--;
                        }
                    }
                }
            }else{ // can not find a largest distance in the remain individual
                break;
            }
        }

        //update t_
        t_ = (double)KneePoint.size()/(double)pop.size();
        return KneePoint;

    }


}
