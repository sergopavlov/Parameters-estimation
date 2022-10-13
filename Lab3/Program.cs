namespace Lab3
{
    internal class Program
    {
        static void Main(string[] args)
        {
            Random random = new();
            double[] coeffstrue = new double[11];
            for (int i = 0; i < 11; i++)
            {
                coeffstrue[i] = i + 1;
            }
            int n = 11;
            double[] points = new double[n];
            for (int i = 0; i < n; i++)
            {
                points[i] = -10 + (20) * random.NextDouble();
                Console.WriteLine(points[i]);
            }
            Func<Speciment, double> f = (s) =>
             {
                 double res = 0;
                 for (int i = 0; i < n; i++)
                 {
                     double p1 = 0;//истинный
                     double p2 = 0;//полученный
                     for (int j = 0; j < 11; j++)
                     {
                         p1 += Math.Pow(points[i], j) * coeffstrue[i];
                         p2 += Math.Pow(points[i], j) * (s.geneList[i] as GeneDouble).value;
                     }
                     res += (p1 - p2) * (p1 - p2);
                 }
                 return res;
             };
            GeneticProblem problem = new GeneticProblem(f, 30000, 11, 0.15, 0.05);
            problem.Solve(10000, 1e-10);

        }
    }

    public abstract class Gene
    {
        public abstract void Mutate(double prob, double val);
    }

    public class GeneDouble : Gene
    {
        static double x0 = 0;
        static double x1 = 20;
        public double value;
        Random random;
        public GeneDouble(double value, Random random)
        {
            this.value = value;
            this.random = random;
        }
        public GeneDouble(Random random)
        {
            this.random = random;
            this.value = x0 + (x1 - x0) * random.NextDouble();
        }
        public override string ToString()
        {
            return value.ToString("0000.00");
        }
        public override void Mutate(double prob, double val)
        {
            if (random.NextDouble() < prob)
            {
                if (random.Next(2) == 0)
                    value *= (1 + val);
                else
                    value *= (1 - val);
            }
        }
    }

    public class Speciment
    {
        public List<Gene> geneList = new List<Gene>();
        int Genecount;
        Random random;
        public void Mutate(double prob, double val)
        {
            foreach (var item in geneList)
            {
                item.Mutate(prob, val);
            }
        }
        public override string ToString()
        {
            System.Text.StringBuilder str = new();
            foreach (var item in geneList)
            {
                str.Append(item.ToString());
                str.Append(" ");
            }
            return str.ToString();
        }
        public Speciment(int genecount, Random random)
        {
            Genecount = genecount;
            this.random = random;
            for (int i = 0; i < genecount; i++)
            {
                geneList.Add(new GeneDouble(random));
            }
        }

        public Speciment(List<Gene> geneList, int genecount, Random random)
        {
            this.geneList = geneList;
            Genecount = genecount;
            this.random = random;
        }

        public static Speciment operator ^(Speciment a, Speciment b)
        {
            int index = a.random.Next(a.Genecount);
            List<Gene> genes = new();
            for (int i = 0; i < index; i++)
            {
                genes.Add(a.geneList[i]);
            }
            for (int i = index; i < a.Genecount; i++)
            {
                genes.Add(b.geneList[i]);
            }
            return new Speciment(genes, a.Genecount, a.random);
        }
    }
    public class GeneticProblem
    {
        Random random = new Random();
        List<Speciment> population = new();
        List<double> Fs = new();
        List<double> probs = new();
        Func<Speciment, double> func;
        int populationSize;
        int genecount;
        double mutatefreq;
        double mutateval;

        public GeneticProblem(Func<Speciment, double> func, int populationSize, int genecount, double mutatefreq, double mutateval)
        {
            this.func = func;
            this.populationSize = populationSize;
            this.genecount = genecount;
            this.mutatefreq = mutatefreq;
            this.mutateval = mutateval;
        }

        void CreateInitialPopulation()
        {
            for (int i = 0; i < populationSize; i++)
            {
                population.Add(new Speciment(genecount, random));
                Fs.Add(0);
                probs.Add(0);
            }
            probs.Add(0);
        }
        int CalcF()
        {
            double summ = 0;
            double summ1 = 0;
            int index = 0;
            for (int i = 0; i < populationSize; i++)
            {
                Fs[i] = func(population[i]);
                probs[i + 1] = 1/Fs[i];
                summ += 1/Fs[i];
                summ1 += Fs[i];
                if (Fs[i] < Fs[index])
                    index = i;
            }
            for (int i = 0; i <= populationSize; i++)
            {
                probs[i] /= summ;
            }
            for (int i = 1; i <= populationSize; i++)
            {
                probs[i] += probs[i - 1];
            }
            Console.WriteLine(summ1 / populationSize);
            return index;
        }
        List<Speciment> CreateNewPopulation()
        {
            List<Speciment> res = new();
            for (int i = 0; i < populationSize; i++)
            {
                double val = random.NextDouble();
                int p1 = 0;
                while (probs[p1 + 1] < val)
                {
                    p1++;
                }
                val = random.NextDouble();
                int p2 = 0;
                while (probs[p2 + 1] < val)
                {
                    p2++;
                }
                res.Add(population[p1] ^ population[p2]);
            }
            return res;
        }
        void Mutate()
        {
            for (int i = 0; i < populationSize; i++)
            {
                population[i].Mutate(mutatefreq, mutateval);
            }
        }
        public void Solve(int maxiter, double eps)
        {
            CreateInitialPopulation();
            int k = 0;
            double minval = 1e+30;
            while (k < maxiter && minval > eps)
            {
                int index = CalcF();
                minval = func(population[index]);
                Console.WriteLine($"{k} {minval} {population[index].ToString()}");
                population = CreateNewPopulation();
                Mutate();
                k++;
            }
        }

    }
}