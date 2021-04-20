package etmo.metaheuristics.moea_ac;

import etmo.core.Solution;
import etmo.core.SolutionSet;

public class MaxSimilarBasedClustering {
	SolutionSet solutionSet_;
	int objNumber_;
	public MaxSimilarBasedClustering(SolutionSet list, int number){
		this.solutionSet_ = list;
		this.objNumber_ = number;
	}
	
	public SolutionSet decomposePopulation(int individualSize){
		int size = solutionSet_.size();
		
		SolutionSet subSet  = new SolutionSet(individualSize);
		
		boolean[] flag =  new boolean[size];
		for(int f=0;f<size;f++){
			flag[f] = false;
		}
		
		int id_max1 = 0;
		int id_max2 = 1;
		double maxDist = computeDistance(solutionSet_.get(0), solutionSet_.get(1));
		for(int k=0; k<size; k++){
			for(int s=k+1;s<size;s++){
				double dis = computeDistance(solutionSet_.get(k), solutionSet_.get(s));
				if(dis > maxDist){
					maxDist = dis;
					id_max1 = k;
					id_max2 = s;
				}
			}
		}//for k
		flag[id_max1] = true;
		flag[id_max2] = true;
		subSet.add(solutionSet_.get(id_max1));
		subSet.add(solutionSet_.get(id_max2));
		
		double[] distances = new double[size];
		int index[] = new int[size];
	 /*compute the distance between each not removed solution in solutionSet_ and subSet[0]*/
		for(int i=0;i<size;i++){
			distances[i] = -1.0;
			Solution sol2 = solutionSet_.get(i);
			if(!flag[i]){
				double minDist = computeDistance(sol2, subSet.get(0));
				int minIndex = 0;
				for(int j=1;j<subSet.size();j++){
					Solution sol3 = subSet.get(j);
					double dist = computeDistance(sol2, sol3);
					if(dist < minDist){
						minDist = dist;
						minIndex = j;
					}
				}//for j
				distances[i] = minDist;
//				System.out.printf("%f\n",minDist);
				if (Double.isNaN(minDist)) System.out.println("distances is NaN!!");
				index[i] = minIndex;
			}//if
		}//for i
		
		//int remain = individualSize - objNumber_;
		int remain = individualSize - 2;
		while(remain > 0){
		 /*find the current solution with the maximum angle to subSets[0]*/
			double maxD = -1.0e+30;
			int maxDisID = -1;
			for(int a=0;a<size;a++){
				if(!flag[a]){
					if(distances[a] > maxD){
						maxD = distances[a];
						maxDisID = a;
					}
				}
			}//for a
		/*maximum angle based addition*/
			flag[maxDisID] = true;
			subSet.add(solutionSet_.get(maxDisID));
		/*update angles*/
		    for(int b=0;b<size;b++){
		    	Solution sol4 = solutionSet_.get(b);
		    	if(!flag[b]){
		    		Solution sol5 = solutionSet_.get(maxDisID);
		    		double angle = computeDistance(sol4,sol5);
		    		if(angle < distances[b]){
		    			distances[b] = angle;
		    			index[b] = subSet.size()-1;
		    		}
		    	}
		    }
		    remain--;
		}//while
		/*SolutionSet[] subSets = new SolutionSet[individualSize];
		for(int m=0;m<individualSize;m++){
			subSets[m] = new SolutionSet();
		}
		
		for(int i=0;i<individualSize;i++){
			subSets[i].add(subSet.get(i));
		}
		for(int i=0;i<size;i++){
			if(!flag[i]){
				double minS = computeDistance(solutionSet_.get(i),subSet.get(0));
			    int mindd = 0;
			    for(int t=1;t<individualSize;t++){
			    	double ss = computeDistance(solutionSet_.get(i),subSet.get(t));
			    	if(ss < minS){
			    		minS = ss;
			    		mindd = t;
			    	}
			    }
				subSets[mindd].add(solutionSet_.get(i));
			}
		}*/
		return subSet;
	}
	
	
	public double computeDistance(Solution so1, Solution so2){
		double dis = 0.0;
		double innerProduc = 0.0;
		
		for(int i=0; i<objNumber_; i++){
			innerProduc += Math.pow(so1.getIthTranslatedObjective(i)- so2.getIthTranslatedObjective(i), 2);
		}
		dis = Math.sqrt(innerProduc);
		if (Double.isNaN(dis)) System.out.println("dis is NaN!!");
		return dis;
	}
}
