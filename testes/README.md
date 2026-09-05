
## Examples

| Directory | Problem type | What it demonstrates |
|-----------|--------------|----------------------|
| [`monodomain`](monodomain/) | Electrophysiology | Monodomain propagation of the action potential. |
| [`eikonal`](eikonal/) | Electrophysiology (fast) | Activation times via the eikonal approximation. |
| [`passive`](passive/) | Mechanics | Passive (unstimulated) tissue mechanics. |
| [`active`](active/) | Mechanics | Active contraction driven by the cell model. |
| [`mechanics`](mechanics/) | Mechanics | Mechanics test cases / benchmarks. |
| [`electromech`](electromech/) | Coupled | Fully coupled electromechanical simulation. |

## Running the Examples

You need to point to the binary dirs, typically in build/app/

### monodomain

Cardiac electrophysiology benchmark

./monodomain -f benchmark_slab_0p5mm.xml -dt 0.001 -t 40.0 -pr 1.0 -c TT2 -m ExplicitEuler
./monodomain -f benchmark_slab_0p2mm.xml -dt 0.001 -t 40.0 -pr 1.0 -c TT2 -m ExplicitEuler
./monodomain -f benchmark_slab_0p1mm.xml -dt 0.001 -t 40.0 -pr 1.0 -c TT2 -m ExplicitEuler

./monodomain -f benchmark_slab_0p1mm.xml -dt 0.0005 -t 60.0 -pr 1.0 -c TT2 -m ExplicitEuler -ksp_type cg -pc_type ilu
./monodomain -f benchmark_slab_0p1mm.xml -dt 0.0005 -t 60.0 -pr 1.0 -c TT2 -m ExplicitEuler -ksp_type cg -pc_type ilu

