import os
import gzip
import argparse

def main():
    parser = argparse.ArgumentParser(description="Split ChromHMM segmentations.bed.gz file into fixed 1000bp windows.")
    parser.add_argument("-i", "--input", required=True, help="Input folder path (including the .segments.bed.gz file).")
    parser.add_argument("-o", "--output", help="Output folder path (split_1k is created in the input directory by default)")
    
    args = parser.parse_args()

    input_folder = args.input
    output_folder = args.output if args.output else os.path.join(input_folder, "split_1k")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print(f"Start processing: {input_folder}")

    for root, dirs, files in os.walk(input_folder):
        if os.path.abspath(root) == os.path.abspath(output_folder):
            continue
            
        for file in files:
            if file.endswith("segments.bed.gz"):
                input_file_path = os.path.join(root, file)
                output_file_name = file.replace(".bed.gz", "_1k.bed")
                output_file_path = os.path.join(output_folder, output_file_name)
                
                # print(f"正在转换: {file} -> {output_file_name}")
                
                with gzip.open(input_file_path, 'rt') as f_in, open(output_file_path, 'w') as f_out:
                    for line in f_in:
                        line = line.strip()
                        if not line: continue
                        
                        parts = line.split('\t')
                        if len(parts) < 4: continue
                        
                        chrom = parts[0]
                        try:
                            start = int(parts[1])
                            end = int(parts[2])
                        except ValueError:
                            continue
                            
                        state = parts[3]
                        length = end - start
                        
                        num_bins = length // 1000
                        if num_bins > 0:
                            for i in range(num_bins):
                                new_start = start + i * 1000
                                new_end = start + (i + 1) * 1000
                                f_out.write(f"{chrom}\t{new_start}\t{new_end}\t{state}\n")
                        # Note: The remaining portion (length % 1000) less than 1000bp is usually discarded according to ChromHMM convention.

    print(f"\nFinish! result saved in: {output_folder}")

if __name__ == "__main__":
    main()