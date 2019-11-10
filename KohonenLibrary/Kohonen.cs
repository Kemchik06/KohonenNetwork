using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace KohonenLibrary
{
    public class Neuron
    {
        private int number;
        public List<double> weights;      
        public double power;
       
        public int Input
        {
            get { return number; }
            set { this.number = value; }
        }

        public List<double> Weights
        {
            get
            { return weights; }

        }

        public double GetRandNum(double min, double max)
        {
            Random random = new Random();
            return random.NextDouble() * (max - min) + min;
        }
        public Neuron(int number, int inputs,double power)
        {
            weights = new List<double>();
            this.number = number;
            this.power = power;
            

            for (int i = 0; i < inputs; i++)
            {

                weights.Insert(i, 0.45);

            }

        }
        
    }
    public class NeuroNet
    {
        private int inputsCount;//колво входов
        public List<double[]> inputs;// список входных данных
        public Neuron[] numsNeur;
        public Neuron this[int index]
        {
            get { return numsNeur[index]; }
            set { numsNeur[index] = value; }
        }
        public List<double[]> Inputs
        {
            get { return inputs; }

        }
        public void SetInputs(double[] input)
        {
            inputs.Add(input);
        }
        public int InputsCount
        {
            get { return inputsCount; }
            set { this.inputsCount = value; }
        }
      
        public NeuroNet(int inputsCount, int neurons_)
        {
            numsNeur = new Neuron[neurons_];
            inputs = new List<double[]>();
            this.inputsCount = inputsCount;
            double power = 1 / neurons_;
            for (int i = 0; i < neurons_; i++)
            {
                Neuron neuron = new Neuron(i, inputsCount,power);

                numsNeur[i] = neuron;
                numsNeur[i].power = power;
            }

        }
        public double EuclideanDistance(Neuron neuron, double[] inputVector)// Евклидова мера между вектором весов и входным вектором
        {
            double Sum = 0;
            for (int i = 0; i < inputVector.Length; i++)
            {
                Sum += (inputVector[i] - neuron.Weights[i]) * (inputVector[i] - neuron.Weights[i]);
            }
            return Math.Sqrt(Sum);
        }
        public void InitPower(int n)
        {
            for(int i = 0; i < numsNeur.Length; i++)
            {
                numsNeur[i].power = 1.0 / n;
            }
        }
        public double[,] Study(int n, int inputsCount)
        {

           
            double minDistance=0;
            int bmuIndex = 0; 
            double[,] newWeights = new double[numsNeur.Length, inputs[0].Length];
           
            double pmin = 0;//в начале потенциал 0 чтобы все могли поучавсвовать в бортьбе

               
                Console.WriteLine(n + " epoha");

          
                if (n > 1 || inputsCount > 0)
                pmin = 0.75;

                    for (int i = 0; i < numsNeur.Length; i++)
                    {
                        if (numsNeur[i].power >= pmin)//нахожу начальное минимальное значение для сравнения с остальными, т е у первого любого попавшегося элемента 
                        {
                            minDistance = EuclideanDistance(this.numsNeur[i], inputs[inputsCount]);//минимальное растояние между нейроном и входом
                            bmuIndex = i;
                            break;
                        }
                    }

                    for (int i = 0; i < numsNeur.Length; i++)
                    {
                        if (numsNeur[i].power >= pmin)
                        {
                            double temp_ED = EuclideanDistance(this.numsNeur[i], inputs[inputsCount]); //находим Евклидово расстояние между i-ым нейроном и k-ым входным вектором 

                            if (temp_ED < minDistance) // если Евклидово расстояние минимально, то это нейрон-победитель 
                            {
                                bmuIndex = i; // индекс нейрона-победителя 
                                minDistance = temp_ED;
                            }
                        }


                    }

                    for (int num = 0; num < this.numsNeur.Length; num++)
                    {
                        if (num == bmuIndex)

                            numsNeur[num].power = numsNeur[num].power - pmin;
                        else
                        {
                            numsNeur[num].power = numsNeur[num].power + (1.0) / numsNeur.Length;
                            if (numsNeur[num].power > (1.0))
                                numsNeur[num].power = 1.0;
                        }
                        if (numsNeur[num].power < (0.0))
                            numsNeur[num].power = 0.0;

                    }
                    for (int i = 0; i < numsNeur.Length; i++)
                    {
                        Console.WriteLine("потенциал " + i + "-го нейрона после " + inputsCount + "-го входа -- " + numsNeur[i].power);
                    }
                    // Console.WriteLine("Веса " + bmuIndex + "-го нейрона"); 
                    
                    for (int j = 0; j < inputs[inputsCount].Length; j++)//изменяем веса нейрона-победителя 
                    {

                        this.numsNeur[bmuIndex].Weights[j] = this.numsNeur[bmuIndex].Weights[j] + LearningConst2(inputsCount) * (inputs[inputsCount][j] - this.numsNeur[bmuIndex].Weights[j]);


                        //Console.WriteLine(j + "-ый вес " + bmuIndex + "-го нейрона: " + this.neurons[bmuIndex].Weights[j]); 
                    }
                 
                
                for (int i = 0; i < numsNeur.Length; i++)//цикл записывающий новые веса
                {
                    for (int j = 0; j < inputs[0].Length; j++)
                    {
                        newWeights[i, j] = numsNeur[i].Weights[j];
                    }
                }
                Console.WriteLine("Веса послей " + n + "-ой эпохи");
                for (int i = 0; i < numsNeur.Length; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        Console.Write(newWeights[i, j] + " , ");
                    }
                    Console.WriteLine();
                }
            
            return newWeights;
                
        }
        public double [,] FinishStudy(int n,int inputsCount)
        {
            double euclDist = 0;
            double[,] newWeights = new double[numsNeur.Length, inputs[0].Length];

                double minDistance = EuclideanDistance(this.numsNeur[0], inputs[inputsCount]);//минимальное растояние между нейроном и входом
                int bmuIndex = 0;//индекс нейрона победителя 
                for (int i = 1; i < this.numsNeur.Length; i++)
                {
                    double temp_ED = EuclideanDistance(this.numsNeur[i], inputs[inputsCount]); //находим Евклидово расстояние между i-ым нейроном и k-ым входным вектором 
                    if (temp_ED < minDistance) // если Евклидово расстояние минимально, то это нейрон-победитель 
                    {
                        bmuIndex = i; // индекс нейрона-победителя 
                        minDistance = temp_ED;
                    }
                }

                // Console.WriteLine("Веса " + bmuIndex + "-го нейрона"); 
                //for (int number = 0; number < numsNeur.Length; number++)
                //{
                    for (int j = 0; j < inputs[inputsCount].Length; j++)//изменяем веса нейрона-победителя 
                    {

                        this.numsNeur[bmuIndex].Weights[j] = this.numsNeur[bmuIndex].Weights[j] +  LearningConst(inputsCount) * (inputs[inputsCount][j] - this.numsNeur[bmuIndex].Weights[j]);

                        //Console.WriteLine(j + "-ый вес " + bmuIndex + "-го нейрона: " + this.neurons[bmuIndex].Weights[j]); 
                    }

                //}
                //euclDist = EuclideanDistance(numsNeur[bmuIndex], inputs[k]);
                //um += euclDist * euclDist;
               


                // error = Math.Sqrt(sum) / inputs.Count;
                // Console.WriteLine(error + " Ошибка , " + inputs.Count);
            //}
            for (int i = 0; i < numsNeur.Length; i++)//цикл записывающий новые веса
            {
                for (int j = 0; j < inputs[0].Length; j++)
                {
                    newWeights[i, j] = numsNeur[i].Weights[j];
                }
            }
            Console.WriteLine("Веса послей " + n + "-ой эпохи");
            for (int i = 0; i < numsNeur.Length; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    Console.Write(newWeights[i, j] + " , ");
                }
                Console.WriteLine();
            }

            euclDist = EuclideanDistance(numsNeur[bmuIndex], inputs[inputsCount]);
            
          
            return newWeights;


        }
        

        public double[,] SecondStudy(int n, int inputsCount)
        {
          
            int time = 1;
            // for (int n = 1; n < 50; n++)
           // {
            double[,] newWeights = new double[numsNeur.Length, inputs[0].Length];

           // Console.WriteLine(n + " epoha");

                ///for (int k = 0; k <inputs.Count; k++) // цикл, в котором предъявляем сети входные вектора - InputVector 
                //{
                    double minDistance = EuclideanDistance(this.numsNeur[0], inputs[inputsCount]);//минимальное растояние между нейроном и входом
                    int bmuIndex = 0;//индекс нейрона победителя 
                    for (int i = 1; i < this.numsNeur.Length; i++)
                    {
                        double temp_ED = EuclideanDistance(this.numsNeur[i], inputs[ inputsCount]); //находим Евклидово расстояние между i-ым нейроном и k-ым входным вектором 
                        if (temp_ED < minDistance) // если Евклидово расстояние минимально, то это нейрон-победитель 
                        {
                            bmuIndex = i; // индекс нейрона-победителя 
                            minDistance = temp_ED;
                        }
                    }

                   // Console.WriteLine("Веса " + bmuIndex + "-го нейрона"); 
                    for (int number = 0; number < numsNeur.Length; number++)
                    {
                        for (int j = 0; j < inputs[inputsCount].Length; j++)//изменяем веса нейрона-победителя 
                        {

                            this.numsNeur[number].Weights[j] = this.numsNeur[number].Weights[j] + NebFunc(n,bmuIndex,number, inputs[inputsCount].Length) *  LearningConst(n) * (inputs[inputsCount][j] - this.numsNeur[number].Weights[j]);

                        //Console.WriteLine(j + "-ый вес " + bmuIndex + "-го нейрона: " + this.neurons[bmuIndex].Weights[j]); 
                        }

                    }
                    //euclDist = EuclideanDistance(numsNeur[bmuIndex], inputs[k]);
                    //sum += euclDist * euclDist;
                    time++;
                

                // error = Math.Sqrt(sum) / inputs.Count;
                // Console.WriteLine(error + " Ошибка , " + inputs.Count);
            //}
            for (int i = 0; i < numsNeur.Length; i++)//цикл записывающий новые веса
            {
                for (int j = 0; j < inputs[0].Length; j++)
                {
                    newWeights[i, j] = numsNeur[i].Weights[j];
                }
            }
            Console.WriteLine("Веса послей " + n + "-ой эпохи");
            for (int i = 0; i < numsNeur.Length; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    Console.Write(newWeights[i, j] + " , ");
                }
                Console.WriteLine();
            }



            return newWeights;
                

           // }
        }

        
        public List<double[]> SwapInput(List<double[]> inputs)
        {
            Random rand = new Random();
            for (int i = inputs.Count - 1; i >= 0; i--)
            {
                int j = rand.Next(i);
                double[] temp = inputs[i];
                inputs[i] = inputs[j];
                inputs[j] = temp;
            }
            return inputs;
        }
        //методы для нормализации даннх
        public double DistOfInput(double[] input)
        {
            double sumOfInput = 0;
            for(int i = 0; i < input.Length;i++)
            {
                sumOfInput += input[i] * input[i];
            }
            sumOfInput = Math.Sqrt(sumOfInput);
            return sumOfInput;
        }
        public double[] NormalizationOfInputs(double[] input)
        {
            double[] norminput = new double [input.Length];
            for(int i = 0; i < input.Length; i++)
            {
                norminput[i] = input[i] / DistOfInput(input);
            }
            return norminput;
        }
        //

            //нахождение нейрона победителя при помощи взвешенной суммы
        public double VzveshSum(Neuron neuron, double[] inputVector)
        {
            double Sum = 0;
            for(int i = 0; i < inputVector.Length; i++)
            {
                Sum += inputVector[i] * neuron.Weights[i];
            }

            return Sum;
        }
        
        public double NebFunc(int t,int bmu,int numOfNuer,int k)
        {
            double dist = Distance(bmu, numOfNuer,k);
            double r = 1/Math.Exp(1/t*t);
            double h = Math.Exp(-(dist / r));
            return h;
        }
        public double Distance(int bmu, int numOfNeur, int k)
        {
            double dist=0;
            for(int j = 0; j < inputs[k].Length; j++)
            {
                dist += (numsNeur[bmu].Weights[j] - numsNeur[numOfNeur].Weights[j]) * (numsNeur[bmu].Weights[j] - numsNeur[numOfNeur].Weights[j]);
            }
            return Math.Sqrt(dist);//расстояние между bmu и тек нейроном
        }
        public double LearningConst(int t)
        {
            double f = 0.1;
            int c = 1000;
            double learningconst = f * Math.Exp(-t / c);
            return learningconst;
        }
        public double LearningConst2(int t)
        {
            double f = 0.5;
            int c = 1000;
            double learningconst = f * Math.Exp(-t / c);
            return learningconst;
        }
        /*public double neiborhooodFunc(int t,int p,int BMU)
        {
            double dist = Distance(BMU, p);
            double r0 = this.numsNeur.Length/2;
            double c = 1000 / Math.Log(r0, 2);
            double r = r0 * Math.Exp((-t) / c);
            double h = Math.Exp((-dist * dist)/2*r*r);
            return h;
        }*/

       

        public double Coordinate (int numofNeur,int numOfCoord)//возращает одну координату опред кластера, для удобства поостроения
        {
            
            return this.numsNeur[numofNeur].Weights[numOfCoord] ;
        }
        
    }
}
