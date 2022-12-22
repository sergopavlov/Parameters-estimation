using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab2v2
{
    public class Parameter
    {
        public string Name { get; init; }
        public double Value { get; set; }
        public override string ToString()
        {
            return $"{Name}: {Value}";
        }
    }
}
