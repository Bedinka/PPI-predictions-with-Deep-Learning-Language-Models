def open_files():
    original_data = open("./fasta_pred.381_250.txt", "r")
    predicted_data = open("./Train18_pred.txt", "r")
    return original_data, predicted_data

def read_files(original_data, predicted_data):
    original_seq = original_data.readlines()
    predicted_seq = predicted_data.readlines()
    return original_seq, predicted_seq

def match_score(original_seq, predicted_seq):
    scores = []
    for original, predicted in zip(original_seq, predicted_seq):
        original = original.strip()
        predicted = predicted.strip()
        score = sum(1 for a, b in zip(original, predicted) if a == b) / len(original) * 100.0
        scores.append(score)
    return scores

if __name__ == "__main__":
    original_data, predicted_data = open_files()
    original_seq, predicted_seq = read_files(original_data, predicted_data)
    original_data.close()
    predicted_data.close()
    scores = match_score(original_seq, predicted_seq)

    with open("./Scores/Train18vsOriginal.txt", "w") as output_file:
        for original, predicted, score in zip(original_seq, predicted_seq, scores):
            output_file.write(f"Original: {original.strip()}\tPredicted: {predicted.strip()}\tScore: {score:.2f}\n")

    print("Comparison scores have been written'.")
