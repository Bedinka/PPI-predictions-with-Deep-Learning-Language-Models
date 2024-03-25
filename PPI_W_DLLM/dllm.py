import data_conversion
import cramming
from transformers import AutoModelForMaskedLM, AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained("pbelcak/UltraFastBERT-1x11-long")
model = AutoModelForMaskedLM.from_pretrained("pbelcak/UltraFastBERT-1x11-long")

vector = data_conversion.main()
encoded_input = tokenizer(vector, return_tensors='pt')
output = model(**encoded_input)