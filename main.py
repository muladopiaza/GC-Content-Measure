from Bio.SeqUtils import gc_fraction
from DNA_CHR01 import dna_chr01
import matplotlib.pyplot as plt
def calculate_windowed_gc(dna_sequence, window_size):
  
  gc_content_list = []
  for i in range(0, len(dna_sequence) - window_size + 1):
    window = dna_sequence[i:i+window_size]
    gc_content = gc_fraction(window) * 100  
    gc_content_list.append(gc_content)
  return gc_content_list


dna_sequence = dna_chr01
window_size = int(len(dna_sequence)/100)

gc_contents = calculate_windowed_gc(dna_sequence, window_size)

print(f"GC content for windows of size {window_size}:")
for i, gc_content in enumerate(gc_contents):
  pass

window_positions = [i for i in range(len(gc_contents))]
plt.plot(window_positions, gc_contents)
plt.xlabel("Window position")
plt.ylabel("GC content (%)")
plt.title(f"GC content variation in {window_size}bp windows")
plt.grid(True)
plt.draw_all(True)
plt.show()
