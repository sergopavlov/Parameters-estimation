namespace Lab3
{
    internal class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello, World!");
        }
    }

    public abstract class Gene
    {

    }

    public class GeneDouble : Gene
    {
        static double x0 = 0;
        static double x1 = 1;
        public double value;
        public GeneDouble(double value)
        {
            this.value = value;
        }
        public GeneDouble(Random random)
        {
            this.value = x0 + (x1 - x0) * random.NextDouble();
        }
    }

    public class Speciment
    {
        private List<Gene> geneList = new List<Gene>();
        int Genecount;
        Random random;

        public Speciment(int genecount, Random random)
        {
            Genecount = genecount;
            this.random = random;
            for (int i = 0; i < genecount; i++)
            {
                geneList.Add(new GeneDouble(random));
            }
        }
    }
    public class GeneticProblem
    {
        Random random;
        List<Speciment> population = new();
        List<double> Fs = new ();
        List<double> probs = new ();
        Func<Speciment, double> func;
        int maxiter;
        double minF;
        int populationSize;
        int genecount;
        void makeIteration()
        {

        }
        void CreateInitialGeneration()
        {
            for (int i = 0; i < populationSize; i++)
            {
                population.Add(new Speciment(genecount, random));
                Fs.Add(0);
                probs.Add(0);
            }
        }
        int CalcF()
        {
            double summ = 0;
            int index = 0;
            for (int i = 0; i < populationSize; i++)
            {
                Fs[i] = func(population[i]);
                probs[i] = Fs[i];
                summ+=Fs[i];
                if (Fs[i] < Fs[index])
                    index = i;
            }
            for (int i = 0; i < populationSize; i++)
            {
                probs[i] /= summ;
            }
            return index;
        }
        List<Speciment> Crossingover()
        {
        }
    }
}