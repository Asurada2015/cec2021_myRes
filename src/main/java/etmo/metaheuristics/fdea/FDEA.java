package etmo.metaheuristics.fdea;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import Jama.Matrix;
import etmo.core.*;
//import etmo.metaheuristics.moea_c.Utils;
import etmo.operators.crossover.EGG;
import etmo.operators.crossover.SBXCrossover;
import etmo.operators.selection.RandomSelection;
import etmo.util.Configuration;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.comparators.ObjectiveComparator;
import etmo.util.ranking.NondominatedRanking;
import etmo.util.ranking.Ranking;
import jmetal.metaheuristics.moead.Utils;

public class FDEA extends Algorithm{
	private int populationSize_;
	
	private SolutionSet population_;
	SolutionSet offspringPopulation_;
	SolutionSet union_;
	
	int generations_;

	Operator crossover_;
	Operator mutation_;
	Operator selection_;
	
	private double[] zideal_; //ideal point
	private double[] znadir_;//Nadir point
	double[][] extremePoints_; // extreme points
	
	int T_;
	int[][] neighborhood_;
	
	double[] pValue;
	double[][] w;
	int t=0;
	
	public FDEA(ProblemSet problem) {
		super(problem);
	} // CAEA_Min_Ideal
    
	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");//写到缓冲区
	        bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	public static void printGD(String path,double[][] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区               
	      for (int i = 0; i < GD.length; i++) {
	    	  for(int j=0;j<GD[i].length;j++){
	    		  bw.write(GD[i][j]+" ");//写到缓冲区
	    	  }
	          bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	
	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int maxGenerations_;

		generations_ = 0;

		maxGenerations_ = ((Integer) this.getInputParameter("maxGenerations")).intValue();
		populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
		mutation_ = operators_.get("mutation");
		crossover_ = operators_.get("crossover");
		selection_ = operators_.get("selection");
		
//		T_ = ((Integer) this.getInputParameter("T")).intValue();
		neighborhood_ = new int[populationSize_][T_];
		
//		int interv;
//		if(problemSet_.get(0).getNumberOfObjectives() == 2){
//			interv = 12;
//		}else if(problemSet_.get(0).getNumberOfObjectives() == 3){
//			interv = 20;
//		}else if(problemSet_.get(0).getNumberOfObjectives() == 5){
//			interv = 24;
//		}else if(problemSet_.get(0).getNumberOfObjectives() == 8){
//			interv = 32;
//		}else if(problemSet_.get(0).getNumberOfObjectives() == 10){
//			interv = 40;
//		}else if(problemSet_.get(0).getNumberOfObjectives() == 15){
//			interv = 60;
//		}else{
//			interv = 30;
//		}
//		pValue = new double[maxGenerations_/interv];
//		w = new double[maxGenerations_/interv][problemSet_.get(0).getNumberOfObjectives()];
		
		initPopulation();// initialize the population;
		
		initIdealPoint();  // initialize the ideal point
		
		initNadirPoint();    // initialize the nadir point
		
		initExtremePoints(); // initialize the extreme points
		int nn = 0;
		while (generations_ < maxGenerations_) {
			
			if(generations_ < maxGenerations_){
				reproduction(generations_, maxGenerations_);
			}else{
				reproduction_Neighbor();
			}
			union_ = ((SolutionSet) population_).union(offspringPopulation_);
			population_.clear();
			SolutionSet[] st = getStSolutionSet(union_,populationSize_);
			double p = 1.0;
			SolutionSet[] subPopulation = null;
			//Autodecomposition
				//estimateIdealPoint(st[0]);
			    updateIdealPoint(st[0]);
			    if(st[0].size() < problemSet_.get(0).getNumberOfObjectives()){
			    	updateNadirPoint(st[1]);
			    }else{
			    	updateNadirPoint(st[0]);
			    }
				
				//estimateNadirPoint(st[1]);
				
				normalizationObjective(st[1]);
				
				/*if((generations_)/100 == 0){
					p = estimation_Curvature(st[0]);
				}*/
				if(generations_ > 0.0*maxGenerations_){
					p = estimation_Curvature(st[0]);
				}
				

				/*if(generations_%interv == 0){
				  pValue[nn] = p; 
				  System.out.println("The current curvature is p = "+p);
				  nn++;
			    }*/
				if(st[1].size() == populationSize_){
					population_ = st[1];
				}else{
					mapping(st[1],p);
				    subPopulation = new MostSimilarBasedSampling(st[1], problemSet_.get(0).getNumberOfObjectives())
							.getIdealPointOrientedPopulation(populationSize_);
				    //Elites Selection to preserve convergence
				    getNextPopulation(subPopulation,generations_,maxGenerations_);
				}
				
			generations_++;
		}//while
		//printGD("FDEA_"+problemSet_.get(0).getNumberOfObjectives()+"Obj_"+problemSet_.get(0).getName()+"_Pvalue.txt",pValue);
		//printGD("FDEA_"+problemSet_.get(0).getNumberOfObjectives()+"Obj_"+problemSet_.get(0).getName()+"_Wvalue.txt",w);
		
		Ranking nodominatedRanking = new NondominatedRanking(population_);
		return nodominatedRanking.getSubfront(0);
		//return population_;
	}//execute
	
	public void initPopulation() throws JMException, ClassNotFoundException {

		population_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problemSet_);

			problemSet_.get(0).evaluate(newSolution);
			problemSet_.get(0).evaluateConstraints(newSolution);

			population_.add(newSolution);
		} // for
	} // initPopulation
	
	public void reproduction(int G, int Gmax) throws JMException{
		offspringPopulation_ = new SolutionSet(populationSize_);
		
		if (crossover_.getClass() == SBXCrossover.class){
			Solution[] parents = new Solution[2];
			for (int i = 0; i < (populationSize_); i++) {
				if (G < Gmax) {
					// obtain parents
					parents = (Solution[]) selection_.execute(population_);

					Solution[] offSpring = (Solution[]) crossover_
							.execute(parents);
					mutation_.execute(offSpring[0]);
					problemSet_.get(0).evaluate(offSpring[0]);
					problemSet_.get(0).evaluateConstraints(offSpring[0]);
					offspringPopulation_.add(offSpring[0]);
				} // if
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

	public void reproduction_Neighbor() throws JMException{
		
		getNeighborhood_Population();
		
		offspringPopulation_ = new SolutionSet(populationSize_);
		Solution[] parents = new Solution[2];
		for (int i = 0; i < (populationSize_); i++) {
			// obtain parents
			int r1 = PseudoRandom.randInt(0, populationSize_-1);
			parents[0] = population_.get(r1);
			double rd = PseudoRandom.randDouble();
			if(rd < 0.2){
				int r2 = PseudoRandom.randInt(0, populationSize_-1);
				while(r1 == r2){
					r2 = PseudoRandom.randInt(0, populationSize_-1);
				}
				parents[1] = population_.get(r2);
			}else{
				int r2 = PseudoRandom.randInt(0, T_-1);
				while(r1 == neighborhood_[r1][r2]){
					r2 = PseudoRandom.randInt(0, T_-1);
				}
				parents[1] = population_.get(neighborhood_[r1][r2]);
			}
			//parents = (Solution[]) selection_.execute(population_);

			Solution[] offSpring = (Solution[]) crossover_
					.execute(parents);
			mutation_.execute(offSpring[0]);
			problemSet_.get(0).evaluate(offSpring[0]);
			problemSet_.get(0).evaluateConstraints(offSpring[0]);
			offspringPopulation_.add(offSpring[0]);
		} // for
	}
	
	public void getNeighborhood_Population(){
		double[] x = new double[populationSize_];
		int[] idx = new int[populationSize_];
		
		for(int i=0;i<populationSize_;i++){
			for(int j=0;j<populationSize_;j++){
				x[j]=computeAngle(population_.get(i),population_.get(j));
				idx[j]=j;
			}
			Utils.minFastSort(x, idx, populationSize_, T_);
			System.arraycopy(idx,0,neighborhood_[i],0,T_);
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
//		System.out.println(generations_ + ": font number = " + front.size());
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
				val = (sol.getObjective(j)-zideal_[j])/(znadir_[j]-zideal_[j]);
				sol.setNormalizedObjective(j, val);
			}//for
		}//for
	}
	
    public double computeAngle(Solution so1, Solution so2){
		double angle = 0.0;
		double distanceToidealPoint1 = so1.getDistanceToIdealPoint();
		double distanceToidealPoint2 = so2.getDistanceToIdealPoint();
		double innerProduc = 0.0;
		for(int i=0; i<problemSet_.get(0).getNumberOfObjectives(); i++){
			innerProduc += so1.getNormalizedObjective(i)*so2.getNormalizedObjective(i);
		}
		angle = innerProduc/(distanceToidealPoint1*distanceToidealPoint2);
		double value = Math.abs(angle);
		if(value < 0.0){
			value = 0.0;
		}else if(value > 1.0){
			value = 1.0;
		}
		angle = Math.acos(value);
		return angle;
	}
    
    public double computeDistance(Solution so1, Solution so2){
		double dis = 0.0;
		double innerProduc = 0.0;
		int objNumber_ = problemSet_.get(0).getNumberOfObjectives();
		for(int i=0; i<objNumber_; i++){
			innerProduc += Math.pow(so1.getIthTranslatedObjective(i)-so2.getIthTranslatedObjective(i), 2);
		}
		dis = Math.sqrt(innerProduc);
		return dis;
	}
    
    public void getNextPopulation(SolutionSet[] subPopulation, int gen, int maxGen){
    	int objNumber = problemSet_.get(0).getNumberOfObjectives();
    	SolutionSet[] ReferenceSet = new LeastSimilarBasedSampling(subPopulation[0],objNumber)
    			.getIdealPointOrientedPopulation(objNumber);
//    	subPopulation[0] -> SB = N, [1] -> SR
		//		ReferenceSet[0] -> Sw  ReferenceSet[1] -> SB - Sw
    	subPopulation[0].clear();
    	
    	SolutionSet subSets[] = new SolutionSet[populationSize_];
    	
    	for(int i=0;i<objNumber;i++){
    		subSets[i] = new SolutionSet();
    		subSets[i].add(ReferenceSet[0].get(i));
    		subPopulation[0].add(ReferenceSet[0].get(i));
    	}
    	
		for(int i=0; i<populationSize_-objNumber;i++){
			  subSets[i+objNumber] = new SolutionSet();
			  subSets[i+objNumber].add(ReferenceSet[1].get(i));
			  subPopulation[0].add(ReferenceSet[1].get(i));
		}

//		subSets -> S1,S2...SN ; subPopulation[0] -> SW U SB = N ; subPopulation[1] -> SR

		// SR -> S1,s2...SN
		for(int i=0; i<subPopulation[1].size();i++){
			Solution s1 = subPopulation[1].get(i);
			double minAngle = computeDistance(s1,subPopulation[0].get(0));
			int minIndex = 0;
			for(int j=1; j<subPopulation[0].size();j++){
				Solution s2 = subPopulation[0].get(j);
				double angle = computeDistance(s1,s2);
				if(angle<minAngle){
				   minAngle = angle;
				   minIndex = j;
				}
			}
			subSets[minIndex].add(s1);
		}


		double[][] centers = new double[populationSize_][problemSet_.get(0).getNumberOfObjectives()];
		for(int i=0;i<populationSize_;i++){
			for(int m=0;m<problemSet_.get(0).getNumberOfObjectives();m++){
				double summ = 0;
				for(int j=0;j<subSets[i].size();j++){
					summ += subSets[i].get(j).getUnitHypersphereObjective(m);
				}
				centers[i][m] = summ/subSets[i].size();
			}
		}
		
		double[] re = new double[problemSet_.get(0).getNumberOfObjectives()];
		for(int i=0;i<problemSet_.get(0).getNumberOfObjectives(); i++){
    		double sum=0;
    		for(int j=0;j<populationSize_; j++){
    			sum += centers[j][i];
    			//sum += subPopulation[0].get(j).getIthTranslatedObjective(i);
    			//sum += subPopulation[0].get(j).getNormalizedObjective(i);
    		}
    		sum = sum/(populationSize_);
    		re[i] = sum;
    		//System.out.print(re[i]+" ");
    		//re[i] = 1.0;
    	}
		/*if(gen%interv == 0){
			w[t] = re;
			t++;
		}*/
		for(int i=0; i<objNumber;i++){
			double rd = PseudoRandom.randDouble();
			if(gen > 0.0*maxGen){
				population_.add(subSets[i].get(0));
			}else{
				double minValue; 
				//minValue = subSets[i].get(0).getDistanceToIdealPoint();
				//minValue = subSets[i].get(0).getSumValue();
				minValue = computeWS(subSets[i].get(0), re);
				int minIndex = 0;
				for(int j=1; j<subSets[i].size();j++){
				   double value; 
				   //value = subSets[i].get(j).getDistanceToIdealPoint();
				   //value = subSets[i].get(j).getSumValue();
				   value = computeWS(subSets[i].get(j), re);
				   if(value < minValue){
					  minValue = value;
					  minIndex = j;
				   }
				}
				population_.add(subSets[i].get(minIndex));
			}
		}
		
		for(int i=objNumber; i<populationSize_;i++){
		    double minValue; 
			//minValue = subSets[i].get(0).getDistanceToIdealPoint();
		    //minValue = subSets[i].get(0).getSumValue();
			minValue = computeWS(subSets[i].get(0), re);
			int minIndex = 0;
			for(int j=1; j<subSets[i].size();j++){
			   double value; 
			   //value = subSets[i].get(j).getDistanceToIdealPoint();
			   //value = subSets[i].get(j).getSumValue();
			   value = computeWS(subSets[i].get(j), re);
			   if(value < minValue){
				  minValue = value;
				  minIndex = j;
			   }
			}
			population_.add(subSets[i].get(minIndex));
		}
    }
    
    public double computeWS(Solution so1, double[] re){
    	double value = 0.0;
    	double sum = 0.0;
    	for(int i=0; i<problemSet_.get(0).getNumberOfObjectives(); i++){
    		sum += re[i];
    	}
    	if(sum == 0){
    		sum = 0.000001;
    	}
    	for(int i=0; i<problemSet_.get(0).getNumberOfObjectives(); i++){
    		re[i] = re[i]/sum;
    		value += so1.getNormalizedObjective(i)*re[i];
    	}
    	return value;
    }
    
    public double computeLWS(Solution s, Solution r){
    	double value=0;
    	for(int i=0;i<problemSet_.get(0).getNumberOfObjectives();i++){
    		value+=s.getNormalizedObjective(i)*r.getNormalizedObjective(i);
    	}
    	return value;
    }
    
    public double estimation_Curvature(SolutionSet solutionSet){
    	SolutionSet solSet = solutionSet;
    	double c = 1.0;
    	int size = solSet.size();
    	int numb = problemSet_.get(0).getNumberOfObjectives();
    	
    	double sum = 0.0;
    	double mean = 0.0;
    	double var = 0.0;
    	
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
  
    	
    	/*strategy*/
    	if(mean >= 0){
    		int T2 = 51;
    		double[] p1 = new double[T2];
    		for(int k=0; k<T2; k++){
        		p1[k] = 1.0 + 0.05*k;
        	}
    		double[] E1 = new double[T2];
    		for(int i=0;i<T2;i++){
    			double ss = 0.0;
    			if(size <= 3){
    				for(int n=0;n<size;n++){
        				Solution sol = solSet.get(n);
        				double sumV = 0.0;
        				for(int m=0;m<numb;m++){
        					sumV += Math.pow(sol.getNormalizedObjective(m), p1[i]);
            			}
        				//sumV = Math.pow(sumV, 1.0/p1[i]);
        				ss += sumV;
        			}
        			E1[i] = (ss)/(size);
    			}else{
    				double minSum = Double.MAX_VALUE;
    				double maxSum = Double.MIN_VALUE;
    				for(int n=0;n<size;n++){
        				Solution sol = solSet.get(n);
        				double sumV = 0.0;
        				for(int m=0;m<numb;m++){
            				sumV += Math.pow(sol.getNormalizedObjective(m), p1[i]);
            			}
        				//sumV = Math.pow(sumV, 1.0/p1[i]);
        				ss += sumV;
        				if(sumV < minSum){
        					minSum = sumV;
        				}
        				if(sumV > maxSum){
        					maxSum = sumV;
        				}
        			}
        			//E1[i] = (ss-minSum-maxSum)/(size-2);
        			E1[i] = (ss)/(size);
    			}
    		}
    		int minID = 0;
        	double min = Math.abs(1.0 - E1[0]);
        	for(int i=1;i<T2;i++){
        		double value1 = Math.abs(1.0 - E1[i]);
        		if(value1 < min){
        			min = value1;
        			minID = i;
        		}
        	}
        	c = p1[minID];
    	}else{
    		int T1 = 17;
    		double[] p2 = new double[T1];
    		for(int k=0; k<T1; k++){
        		p2[k] = 1.0 - 0.03*(k);
        	}
    		double[] E2 = new double[T1];
    		for(int i=0;i<T1;i++){
    			double ss = 0.0;
    			if(size <= 3){
    				for(int n=0;n<size;n++){
        				Solution sol = solSet.get(n);
        				double sumV = 0.0;
        				for(int m=0;m<numb;m++){
        					sumV += Math.pow(sol.getNormalizedObjective(m), p2[i]);
            			}
        				sumV = Math.pow(sumV, 1.0/p2[i]);
        				ss += sumV;
        			}
        			E2[i] = (ss)/(size);
    			}else{
    				double minSum = Double.MAX_VALUE;
    				double maxSum = Double.MIN_VALUE;
    				for(int n=0;n<size;n++){
        				Solution sol = solSet.get(n);
        				double sumV = 0.0;
        				for(int m=0;m<numb;m++){
            				sumV += Math.pow(sol.getNormalizedObjective(m), p2[i]);
            			}
        				sumV = Math.pow(sumV, 1.0/p2[i]);
        				ss += sumV;
        				if(sumV < minSum){
        					minSum = sumV;
        				}
        				if(sumV > maxSum){
        					maxSum = sumV;
        				}
        			}
        			//E2[i] = (ss-minSum-maxSum)/(size-2);
        			E2[i] = (ss)/(size);
    			}
    		}
    		int minID = 0;
        	double min = Math.abs(1.0 - E2[0]);
        	for(int i=1;i<T1;i++){
        		double value1 = Math.abs(1.0 - E2[i]);
        		if(value1 < min){
        			min = value1;
        			minID = i;
        		}
        	}
        	c = p2[minID];
        	double c1, c2;
    	}

    	//c = PseudoRandom.randDouble(p[minID[1]],p[minID[0]]);
   
    	if(c!=1.0){
        	if(cv < 0.15){
        		double rd = PseudoRandom.randDouble();
        		if(mean < 0 && rd < 0.95){
        			c = 1.0 - cv;
        		}else if(mean > 0 && rd < 0.95){
        			c = 1.0 + cv;
        		}
        	}/*else if(cv<0.2 && cv>=0.1){
        		double rd = PseudoRandom.randDouble();
        		if(mean < 0 && rd < 0.8){
        			c = 0.88;
        		}else if(mean > 0 && rd < 0.8){
        			c = 1.25;
        		}
        	}*/
    	}
    	//c=1.7;
    	return c;
    }
    
    public void mapping(SolutionSet solSet, double curvature){
    	double p = curvature;
    	int size = solSet.size();
    	int numb = problemSet_.get(0).getNumberOfObjectives();
    	
    	for(int i=0; i<size; i++){
    		Solution sol = solSet.get(i);
    		double normDistance = 0.0;
    		double sumValue = 0.0;
    		double distance = 0.0;
    		for(int j=0; j<numb; j++){
    			normDistance += Math.pow(sol.getNormalizedObjective(j), p);
    			sumValue +=  sol.getNormalizedObjective(j);
    			distance += Math.pow(sol.getNormalizedObjective(j), 2);
    		}
    		normDistance = Math.pow(normDistance, 1/p);
    		distance = Math.sqrt(distance);
    		
    		sol.setDistanceToIdealPoint(normDistance);
    		sol.setDistanceToNadirPoint(distance);
    		sol.setSumValue(sumValue);
    		
    		if(sol.getDistanceToIdealPoint() == 0){
				//System.out.println("Error: This solution is in the origin");
				double dis = 0.0;
				/*for(int j=0; j<numb; j++){
					sol.setNormalizedObjective(j, 0.01);
					dis += Math.pow(sol.getNormalizedObjective(j), p);
				}
				dis = Math.pow(dis, 1/p);*/
				sol.setDistanceToIdealPoint(0.001);
				sol.setDistanceToNadirPoint(0.001);
				//System.exit(0);
			}
    		
    		for(int j=0; j<problemSet_.get(0).getNumberOfObjectives(); j++){
    			double value1 = sol.getNormalizedObjective(j)/sol.getDistanceToIdealPoint();
    			sol.setIthTranslatedObjective(j, value1);
    			double value2 = sol.getNormalizedObjective(j)/sol.getDistanceToNadirPoint();
    			sol.setUnitHypersphereObjective(j, value2);
    		}
    	}
    }
    
    public int permutation(int N, int M){
    	int result = 1;
    	int n = N;
    	int m = M;
    	
    	for(int i=m;i>0;i--){
    		result *= n;
    		n--;
    	}
    	
    	return result;
    }
    public int combination(int N, int M){
    	int result = 1;
    	int n = N;
    	int m = M;
    	
    	int half = n/2;
    	if(m > half){
    		m = n - m;
    	}
    	
    	int  numerator = permutation(n, m);
    	
    	int denominator = permutation(m, m);
    	
    	result = numerator/denominator;
    	
    	return result;
    }
    
    public void sss(){
    	int numb = 12;
    	int sumId = 0;
    	double p = 1.0;
    	double sumV = 0.0;
    	for(int m=1;m<=numb;m++){
    		int c = combination(numb,m);
    		sumId += c;
    		sumV += (((Math.pow(m,1.0-(1.0/p))-1.0)/Math.sqrt(numb))*(double)c);
    	}
    	
    	double E = sumV/sumId;
    }

}
