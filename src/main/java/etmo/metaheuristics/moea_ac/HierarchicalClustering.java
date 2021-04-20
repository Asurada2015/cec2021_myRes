package etmo.metaheuristics.moea_ac;

import etmo.core.Solution;
import etmo.core.SolutionSet;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/*
 * 采用非递归的方法来实现层次聚类，大大的降低了时间复杂度
 */
public class HierarchicalClustering {
	List<SolutionSet> list = new <SolutionSet>ArrayList();
	public HierarchicalClustering(List list){
		this.list = list;
	}
	public List<SolutionSet> clusteringAnalysis(int clusteringSize, double p){
		int size = list.size();
		double minAngle = Double.MAX_VALUE;
		double angle = 0.0;
		double[] minAngles = new double[size];//每个类与之最近类之间的角度
		int[] minIndexs = new int[size];//每个类与之最近的类
		int[] index = new int[4];//当前最相似的两个类
		index[0] = index[1] = index[2] = index[3] = -1;
		for(int i=0; i<size; i++){
			list.get(i).setCurve(p);
			double min = Double.MAX_VALUE;
			for(int j=0;j<size;j++){
				if(i != j){
					angle = computeDistance(list.get(i).getAdaptiveCentroid(),list.get(j).getAdaptiveCentroid());
					if (angle ==  Double.NaN) System.out.println("angle is NaN");
					//angle = computeAngle(list.get(i).getCentroid(),list.get(j).getCentroid());
					if(min > angle){
						min = angle;
						index[0] = i;
						index[1] = j;
					}
				}
			}
			minAngles[i] =  min;
			minIndexs[i] = index[1];
			if(minAngle > min){
				minAngle = min;
				index[2] = index[0];
				index[3] = index[1];
			}
		}
		while(size > clusteringSize){
			if (index[2] == -1 && index[3] == -1){
				System.out.println("index[2] == -1 && index[3] == -1");
				System.exit(0);
			}
			SolutionSet sols = (list.get(index[2]).union(list.get(index[3])));
			sols.setCurve(p);
			list.get(index[2]).setRemove(true);	
			list.remove(index[3]);
			list.add(index[3], sols);
			 /*
		     * 更新把index[2]个体当作角度最近个体i的最近个体indexs[i];
		     */
			for(int i=0;i<list.size();i++){
				if(minIndexs[i]==index[2] && !list.get(i).isRemove()){
					double min = Double.MAX_VALUE;
		    		int sb = -1;
		    		double ss = 0.0;
		    		for(int j=0;j<list.size();j++){
		    			if(!list.get(j).isRemove() && i!=j){
		    				//ss = computeAngle(list.get(i).getCentroid(),list.get(j).getCentroid());
		    				ss = computeDistance(list.get(i).getAdaptiveCentroid(),list.get(j).getAdaptiveCentroid());
		    				if(min > ss){
			    				min = ss;
			    				sb = j;
			    			}//if
		    			}//if
		    		}//for
		    		minAngles[i] = min;
			    	minIndexs[i] = sb;
				}//if
			}//for
			 /*
		     * 更新把index[3]个体当作角度最近个体i的最近个体indexs[i];
		     */
			for(int i=0;i<list.size();i++){
				if(minIndexs[i]==index[3] && !list.get(i).isRemove()){
					double min = Double.MAX_VALUE;
		    		int sb = -1;
		    		double ss = 0.0;
		    		for(int j=0;j<list.size();j++){
		    			if(!list.get(j).isRemove() && i!=j){
		    				//ss = computeAngle(list.get(i).getCentroid(),list.get(j).getCentroid());
		    				ss = computeDistance(list.get(i).getAdaptiveCentroid(),list.get(j).getAdaptiveCentroid());
		    				if(min > ss){
			    				min = ss;
			    				sb = j;
			    			}//if
		    			}//if
		    		}//for
		    		minAngles[i] = min;
			    	minIndexs[i] = sb;
				}//if
			}//for
			/*
		     * 更新当前最近角度的两个个体的index;
		     */
			double sAngle = Double.MAX_VALUE;
			for(int k=0;k<list.size();k++){
				if(!list.get(k).isRemove()){
					if(sAngle > minAngles[k]){
						sAngle = minAngles[k];
						index[2] = k;
						index[3] = minIndexs[k]; 	
					}
				}
			}
	
			size--;
		}//while
		/*int sss = 0;
		for(int i=0;i<list.size();i++){
			if(!list.get(i).isRemove()){
				sss++;
			}
		}*/
		/*if(sss != 16){
			System.out.println("sss = " + sss);
			System.exit(0);
		}*/
		Iterator<SolutionSet> iterator = list.iterator();
		while(iterator.hasNext()){
			if(iterator.next().isRemove()){
				iterator.remove();
			}
		}
		return this.list;
	}
	
	/*
     * 求两个个体之间的角度值
     */
	public double computeDistance(Solution s1, Solution s2){
		double distance = 0.0;
		for(int i=0; i<list.get(0).get(0).getNumberOfObjectives(); i++){
			distance += Math.pow((s1.getIthTranslatedObjective(i)-s2.getIthTranslatedObjective(i)), 2);
		}
		distance = Math.sqrt(distance);
		return distance;
	}//computeDistance

}
