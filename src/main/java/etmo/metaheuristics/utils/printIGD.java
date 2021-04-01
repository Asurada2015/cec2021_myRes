package etmo.metaheuristics.utils;

import etmo.util.Configuration;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

public class printIGD {
    public static void printIGDtoText(String path, double[][] cpIGD, int taskNumber, int times){
        try{
            FileOutputStream fos = new FileOutputStream(path, true);
            OutputStreamWriter osw = new OutputStreamWriter(fos);
            BufferedWriter bw = new BufferedWriter(osw);

            for (int i = 0; i < taskNumber; i++) {
                for (int j = 0; j < times; j++)
                    bw.write(Double.toString(cpIGD[i][j]) + " ");
                bw.newLine();
            }
            /* Close the file */
            bw.close();

        } catch (IOException e){
            Configuration.logger_.severe("Error acceding to the file");
            e.printStackTrace();
        }
    }
}
