using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab2v2
{
    public abstract class DirectProblem
    {
        public int ComponentCount;
        public abstract double CalcValueAtPoint(int index, List<Parameter> Params);
        public abstract List<List<double>> DiffComponentByParameter(List<Parameter> Params);

    }
}
