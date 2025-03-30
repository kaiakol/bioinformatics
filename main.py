from processing import preprocess_input, calculate_background_model, run_meme, postprocess_output

#Run software from this file
def main():
    input_file = "crp0.fna"
    output_dir = "meme_output"
    motif_length = (6, 50)
    background_model = "background_model.txt"
    k_range = (2, 10)

    sequences = preprocess_input(input_file)
    calculate_background_model(input_file, background_model)
    run_meme(input_file, output_dir, motif_length, background_model)
    postprocess_output(output_dir, k_range)

if __name__ == "__main__":
    main()
