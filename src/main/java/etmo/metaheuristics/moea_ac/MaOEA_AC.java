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
public class MaOEA_AC extends Algorithm{
	private SolutionSet population_;
	private SolutionSet offspringPopulation_;
	private SolutionSet union_;
	
	private int populationSize_;
	int generations_;
	int maxGenerations_;
	
	Operator selection_;
	Operator crossover_;
	Operator mutation_;
	Operator learning_;
	
	private double[] zideal_; //ideal point
	private double[] znadir_;//Nadir point
	double[][] extremePoints_; // extreme points
	private double p = 1.0;
	
	public MaOEA_AC(ProblemSet problem) {
		super(problem);
		//zideal_ = new double[problem.getNumberOfObjectives()];
		//znadir_ = new double[problem.getNumberOfObjectives()];
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
		 * Enter the main loop£¬into the process of evolution
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
		
		Ranking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);
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
		population_ = new SolutionSet(populationSize_);
		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problemSet_);
			problemSet_.get(0).evaluate(newSolution);
			problemSet_.get(0).evaluateConstraints(newSolution);
			population_.add(newSolution);
		} // for
	} // initPopulation

	/*
	 * step3 and step4:Mating Selection and Recombination to generate Offspring Populations
	 */
	public void generateOffspringPopulation() throws JMException{
		offspringPopulation_ = new SolutionSet(populationSize_);
		if (crossover_.getClass() == SBXCrossover.class){
			Solution[] parents = new Solution[2];
			for (int i = 0; i < (populationSize_); i++) {
				// obtain parents
				parents = (Solution[]) selection_.execute(population_);

				Solution[] offSpring = (Solution[]) crossover_
						.execute(parents);
				mutation_.execute(offSpring[0]);
				problemSet_.get(0).evaluate(offSpring[0]);
				problemSet_.get(0).evaluateConstraints(offSpring[0]);
				offspringPopulation_.add(offSpring[0]);
			} // for
		}
		else if (crossover_.getClass() == EGG.class){
			int[] permutation = new int[populationSize_];
			Utils.randomPermutation(permutation, populationSize_);
			etmo.util.Ranking ranking = new etmo.util.Ranking(population_);
			SolutionSet front = ranking.getSubfront(0);

			front.Suppress();
			SolutionSet KP = null;
			if(front.size()> problemSet_.get(0).getNumberOfObjectives())
				KP = findingKneePoint(front);
			if(KP == null){
				KP = population_;
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
				parents[0] = population_.get(n);
				child = (Solution) crossover_.execute(parents);

				mutation_.execute(child);

				problemSet_.get(0).evaluate(child);
				problemSet_.get(0).evaluateConstraints(child);
				offspringPopulation_.add(child);
			} // for
		}

	}



	/*
	 * step5:Environmental Selection
	 */
	public void environmentalSelection(int G){
		/*
		 * step5.1:Combine the Population and the Offspring Population
		 */
		union_ = ((SolutionSet) population_).union(offspringPopulation_);
		/*
		 * step5.2:Normalization the Combined Population
		 */
		SolutionSet[] st = getStSolutionSet(union_,populationSize_);

//		St  = st[1]
//		ND  = st[0]

		//estimateIdealPoint(st[1]);
		//estimateNadirPoint(st[1]);
		updateIdealPoint(st[0]);
		updateNadirPoint(st[0]);
		normalizationObjective(st[1]);
		if(G >= 0.0*maxGenerations_){
			//p = estimation_Curvature(st[1]);
			p = estimation_Curvature(st[0]);
			//p = 0.5;
		}
		if(G % 50 == 0){
//			System.out.println("p = "+p);
		}
		computeDistanceToIdealPoint(st[1],p);
		population_.clear();


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
		
		bestSolutionSelection(list);	
	}
	
	 /*
		 * Estimate the Ideal Point 
		 */
		public void estimateIdealPoint(SolutionSet solutionSet){
			for(int i=0; i<problemSet_.get(0).getNumberOfObjectives();i++){
				zideal_[i] = 1.0e+30;
				for(int j=0; j<solutionSet.size();j++){
					if(solutionSet.get(j).getObjective(i) < zideal_[i]){
						zideal_[i] = solutionSet.get(j).getObjective(i);
					}//if
				}//for
			}//for
		}
		
		/*
		 * Estimate the Nadir Point 
		 */
	    public void estimateNadirPoint(SolutionSet solutionSet){
	    	for(int i=0; i<problemSet_.get(0).getNumberOfObjectives();i++){
				znadir_[i] = -1.0e+30;
				for(int j=0; j<solutionSet.size();j++){
					if(solutionSet.get(j).getObjective(i) > znadir_[i]){
						znadir_[i] = solutionSet.get(j).getObjective(i);
					}//if
				}//for
			}//for
		}
		
	    void copyObjectiveValues(double[] array, Solution individual) {
			for (int i = 0; i < individual.getNumberOfObjectives(); i++) {
				array[i] = individual.getObjective(i);
			}
		}
		double asfFunction(Solution sol, int j) {
			double max = Double.MIN_VALUE;
			double epsilon = 1.0E-6;

			int obj = problemSet_.get(0).getNumberOfObjectives();

			for (int i = 0; i < obj; i++) {

				double val = Math.abs((sol.getObjective(i) - zideal_[i])
						/ (znadir_[i] - zideal_[i]));

				if (j != i)
					val = val / epsilon;

				if (val > max)
					max = val;
			}

			return max;
		}

		double asfFunction(double[] ref, int j) {
			double max = Double.MIN_VALUE;
			double epsilon = 1.0E-6;

			int obj = problemSet_.get(0).getNumberOfObjectives();

			for (int i = 0; i < obj; i++) {

				double val = Math.abs((ref[i] - zideal_[i])
						/ (znadir_[i] - zideal_[i]));
				

				if (j != i)
					val = val / epsilon;

				if (val > max)
					max = val;
			}

			return max;
		}

		void initIdealPoint() {
			int obj = problemSet_.get(0).getNumberOfObjectives();
			zideal_ = new double[obj];
			for (int j = 0; j < obj; j++) {
				zideal_[j] = Double.MAX_VALUE;

				for (int i = 0; i < population_.size(); i++) {
					if (population_.get(i).getObjective(j) < zideal_[j])
						zideal_[j] = population_.get(i).getObjective(j);
				}
			}
		}

		void updateIdealPoint(SolutionSet pop){
			for (int j = 0; j < problemSet_.get(0).getNumberOfObjectives(); j++) {
				for (int i = 0; i < pop.size(); i++) {
					if (pop.get(i).getObjective(j) < zideal_[j])
						zideal_[j] = pop.get(i).getObjective(j);
				}
			}
		}
		
		void initNadirPoint() {
			int obj = problemSet_.get(0).getNumberOfObjectives();
			znadir_ = new double[obj];
			for (int j = 0; j < obj; j++) {
				znadir_[j] = Double.MIN_VALUE;

				for (int i = 0; i < population_.size(); i++) {
					if (population_.get(i).getObjective(j) > znadir_[j])
						znadir_[j] = population_.get(i).getObjective(j);
				}
			}
		}

		public void initExtremePoints() {
			int obj = problemSet_.get(0).getNumberOfObjectives();
			extremePoints_ = new double[obj][obj];
			for (int i = 0; i < obj; i++){
				for (int j = 0; j < obj; j++){
					extremePoints_[i][j] = 1.0e+30;
				}
			}
			
		}

		void updateNadirPoint(SolutionSet pop){
			
			updateExtremePoints(pop);

			
			int obj = problemSet_.get(0).getNumberOfObjectives();
			double[][] temp = new double[obj][obj];

			for (int i = 0; i < obj; i++) {
				for (int j = 0; j < obj; j++) {
					double val = extremePoints_[i][j] - zideal_[j];
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

					double aj = 1.0 / AL.get(j, 0) + zideal_[j];
			

					if ((aj > zideal_[j]) && (!Double.isInfinite(aj)) && (!Double.isNaN(aj)))
						znadir_[j] = aj;
					else {
						sucess = false;
						break;
					}
				}
			} 
			else 
				sucess = false;
			
			
			if (!sucess){
				double zmax[] = computeMaxPoint(pop);
				for (int j = 0; j < obj; j++) {
					znadir_[j] = zmax[j];
				}
			}
		}

		public void updateExtremePoints(SolutionSet pop){
			for (int i = 0; i < pop.size(); i++)
				updateExtremePoints(pop.get(i));
		}

		public void updateExtremePoints(Solution individual){
			int obj = problemSet_.get(0).getNumberOfObjectives();
			for (int i = 0; i < obj; i++){
				double asf1 = asfFunction(individual, i);
				double asf2 = asfFunction(extremePoints_[i], i);
				
				if (asf1 < asf2){
					copyObjectiveValues(extremePoints_[i], individual);
				}
			}
		}

		double[] computeMaxPoint(SolutionSet pop){
			int obj = problemSet_.get(0).getNumberOfObjectives();
			double zmax[] = new double[obj];
			for (int j = 0; j < obj; j++) {
				zmax[j] = Double.MIN_VALUE;

				for (int i = 0; i < pop.size(); i++) {
					if (pop.get(i).getObjective(j) > zmax[j])
						zmax[j] = pop.get(i).getObjective(j);
				}
			}
			return zmax;
		}
		
	    /*
	     * Normalization
	     */
		public void normalizationObjective(SolutionSet solutionSet){
			for(int i=0; i<solutionSet.size(); i++){
				Solution sol = solutionSet.get(i);
				for(int j=0; j<problemSet_.get(0).getNumberOfObjectives(); j++){
					double val = 0.0;
					val = (sol.getObjective(j)-zideal_[j])/(znadir_[j]-zideal_[j]+1e-30);
					if (val == Double.POSITIVE_INFINITY) System.out.println("set normalization is infinity!!");
					sol.setNormalizedObjective(j, val);
				}//for
			}//for
		}
	
	 /*
     * Compute the Convergence Distance of each Solutions Which use the distance of 
     * each solution to the Ideal Point
     */
    public void computeDistanceToIdealPoint(SolutionSet solutionSet, double p){
    	for(int i=0; i<solutionSet.size(); i++){
    		Solution sol = solutionSet.get(i);
    		double normDistance = 0.0;
    		double sumValue = 0.0;
    		double pValue = 0.0;
    		for(int j=0; j<problemSet_.get(0).getNumberOfObjectives(); j++){
    			normDistance += sol.getNormalizedObjective(j)*sol.getNormalizedObjective(j);
    			sumValue +=  sol.getNormalizedObjective(j);
    			pValue += Math.pow(sol.getNormalizedObjective(j), p);
    		}
    		normDistance = Math.sqrt(normDistance);
    		pValue = Math.pow(pValue, 1.0/p);
    		sol.setDistanceToIdealPoint(normDistance);
//    		debug
    		if (Double.isInfinite(sumValue)) System.out.println("sumValue = INFINITY");
    		sol.setSumValue(sumValue);
    		for(int j=0; j<problemSet_.get(0).getNumberOfObjectives(); j++){
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
    
    public void bestSolutionSelection(List<SolutionSet> list){
    	boolean[] flag = new boolean[list.size()];
    	SolutionSet solSet = new SolutionSet();
       for(int i=0;i<list.size();i++){
    	   for(int j=0;j<list.get(i).size();j++){
    		   list.get(i).get(j).setClusterID(i);
    		   solSet.add(list.get(i).get(j));
    	   }
       }
       SolutionSet cornerSet = new MaxSimilarBasedClustering(solSet, problemSet_.get(0).getNumberOfObjectives())
    		   .decomposePopulation(problemSet_.get(0).getNumberOfObjectives());
       for(int m=0;m<problemSet_.get(0).getNumberOfObjectives();m++){
    	   if(!flag[cornerSet.get(m).getClusterID()]){
    		   flag[cornerSet.get(m).getClusterID()] = true;
    		   population_.add(cornerSet.get(m));
    	   }/*else{
    		   System.out.println("ÖØ¸´Ìí¼ÓÁË");
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
   			population_.add(sols.get(minFitnessID));
    	   }
       }
   }
   
   public double estimation_Curvature(SolutionSet solutionSet){
   	SolutionSet solSet = solutionSet;
   	double c = 1.0;
   	int size = solSet.size();
   	int numb = problemSet_.get(0).getNumberOfObjectives();
   	
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
   	
   	double alpha = 0.5 - problemSet_.get(0).getNumberOfObjectives()*0.02;
   	
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


private SolutionSet findingKneePoint(SolutionSet pop) {
	Matrix Pop_M, ones_M, ExtremePoint_M, Hyperplane, Distance_M;
	SolutionSet KneePoint = new SolutionSet();
	SolutionSet ExtremePoint = new SolutionSet(problemSet_.get(0).getNumberOfObjectives());
	SolutionSet TempPop = new SolutionSet(pop.size());
	double [] f_max = new double [problemSet_.get(0).getNumberOfObjectives()];
	double [] f_min = new double [problemSet_.get(0).getNumberOfObjectives()];
	for (int i=0;i<problemSet_.get(0).getNumberOfObjectives();i++) {
		f_max[i]=pop.get(0).getObjective(i);f_min[i]=pop.get(0).getObjective(i);
	}
	for(int i=0;i<pop.size();i++) {
		TempPop.add(pop.get(i));
		for (int j=0;j<problemSet_.get(0).getNumberOfObjectives();j++) {
			if(pop.get(i).getObjective(j)>f_max[j]) f_max[j]=pop.get(i).getObjective(j);
			if(pop.get(i).getObjective(j)<f_min[j]) f_min[j]=pop.get(i).getObjective(j);
		}
	}
	ones_M = new Matrix(problemSet_.get(0).getNumberOfObjectives(),1,1.0);
	//1. find extrem solution
	for(int i=0;i<problemSet_.get(0).getNumberOfObjectives();i++){
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
	r_ = r_*Math.exp(-((1.0-t_/0.5)/(double)problemSet_.get(0).getNumberOfObjectives()));
	double[] R = new double[problemSet_.get(0).getNumberOfObjectives()];
	for(int i=0;i<problemSet_.get(0).getNumberOfObjectives();i++)
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
					for(int j=0;j<problemSet_.get(0).getNumberOfObjectives();j++){
						if(Math.abs(s.getObjective(j)-currentKneeP.getObjective(j))>R[j]) {flag = false;break;}
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
