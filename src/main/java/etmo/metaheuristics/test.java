package etmo.metaheuristics;

//import Jama.Matrix;

import java.util.Arrays;

public class test {
    public static void main(String[] args) {
        int[] answers = {0,2,0,2,1};
        int n = answers.length;
        Arrays.sort(answers);
        int sum = 0;
        int tmp = 0;
        for (int i = 1; i < n; i++){
            if (answers[i] != answers[i - 1]){
                int len = i - tmp, num = answers[i - 1] + 1;
                sum += (len / num + len % num) * num;
                tmp = i;
            }
            if (i == n - 1){
                int len = i - tmp + 1, num = answers[i] + 1;
                sum += (len / num + len % num) * num;
                tmp = i;
            }
        }


    }
}



