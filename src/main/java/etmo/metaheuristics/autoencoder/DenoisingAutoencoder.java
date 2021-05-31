package etmo.metaheuristics.autoencoder;

import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;

import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.nd4j.linalg.activations.Activation;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.learning.config.AdaGrad;
import org.nd4j.linalg.lossfunctions.LossFunctions;

public class DenoisingAutoencoder {

    private static final int hiddenNum = 20;

    int inN_;
    int outN_;
    int hN_;

    int dN_;

    NDArray input_;
    NDArray output_;

    MultiLayerNetwork net;

    public DenoisingAutoencoder(double[][] features, double[][] labels){
        dN_ = features.length;
        inN_ = features[0].length;
        outN_ = labels[0].length;
        input_ = new NDArray(features);
        output_ = new NDArray(labels);
        hN_ = hiddenNum;
    }

    public DenoisingAutoencoder(NDArray features, NDArray labels){
        dN_ = features.rows();
        inN_ = features.columns();
        outN_ = labels.columns();
        input_ = features;
        output_ = features;
        hN_ = hiddenNum;
    }

    public void train(){
        // 网络初始化
        var conf = new NeuralNetConfiguration.Builder()
                .seed(12345)
                .weightInit(WeightInit.XAVIER)
                .updater(new AdaGrad(0.05))
                .activation(Activation.SIGMOID)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .l2(0.0001)
                .list()
                .layer(0, new DenseLayer.Builder().nIn(inN_).nOut(hN_).build())
                .layer(1, new OutputLayer.Builder().nIn(hN_).nOut(outN_)
                        .lossFunction(LossFunctions.LossFunction.MSE)
                        .build())
                .validateOutputLayerConfig(false)
                .build();

        net = new MultiLayerNetwork(conf);
//        net.setListeners(new ScoreIterationListener(1));

        net.fit(input_, output_);
    }

    public double[][] predict(double[][] input){
        NDArray in = new NDArray(input);
        INDArray out = net.output(in);
        return out.toDoubleMatrix();
    }



    public double[][] predict(NDArray input){
        INDArray out = net.output(input);
        return out.toDoubleMatrix();
    }


//    double[] inputLayer, outputLayer, hiddenLayer;
//    double[][] i2h, h2o;
//    double[] hB, oB;
//    double learningRate;
//
//    public void init(int hiddenNum, double learningRate, Variable[] source, Variable[] target) throws JMException {
//
//        int inputNum = source.length;
//        int outputNum = target.length;
//        // source work for target
//
//        this.inputLayer = new double[outputNum];
//        this.outputLayer = new double[outputNum];
//        this.hiddenLayer = new double[hiddenNum];
//        this.learningRate = learningRate;
//        this.i2h = new double[outputNum][hiddenNum];
//        this.h2o = new double[hiddenNum][outputNum];
//        this.hB = new double[hiddenNum];
//        this.oB = new double[outputNum];
//
//        if (inputNum > outputNum){
//            for (int i = 0; i < outputNum; i++){
//                Arrays.fill(i2h[i], Math.random());
//                inputLayer[i] = source[i].getValue();
//            }
//            for (int i = 0; i < outputNum; i++){
//                Arrays.fill(h2o[i], Math.random());
//                outputLayer[i] = target[i].getValue();
//            }
//        }
//        else{
//            for (int i = 0; i < outputNum; i++){
//                Arrays.fill(i2h[i], Math.random());
//                if (i < inputNum) inputLayer[i] = source[i].getValue();
//                else inputLayer[i] = target[i].getValue();
//            }
//            for (int i = 0; i < outputNum; i++){
//                Arrays.fill(h2o[i], Math.random());
//                outputLayer[i] = target[i].getValue();
//            }
//        }






}






