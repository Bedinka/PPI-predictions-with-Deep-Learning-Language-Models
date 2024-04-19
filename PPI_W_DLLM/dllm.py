import YangLabIntern.PPI_W_DLLM.data_autoencoding as data_autoencoding
import cramming
from transformers import AutoModelForMaskedLM, AutoTokenizer
import numpy as np

tokenizer = AutoTokenizer.from_pretrained("pbelcak/UltraFastBERT-1x11-long")
model = AutoModelForMaskedLM.from_pretrained("pbelcak/UltraFastBERT-1x11-long")

vector  = data_autoencoding.main()
res = np.char.mod('%s', matrix).tolist()
encoded_input = tokenizer(str(vector), return_tensors='pt')
output = model(**encoded_input)
print(output)
#ValueError: text input must be of type `str` (single example), `List[str]` (batch or single pretokenized example) or `List[List[str]]` (batch of pretokenized examples).
