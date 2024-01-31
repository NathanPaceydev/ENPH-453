import pandas as pd
import matplotlib.pyplot as plt
import sys

def process_cnf(file_path):
    # Read the file
    df = pd.read_csv(file_path, sep='\t', skiprows=35, names=['Channel', 'Energy', 'Counts', 'Rate'])

    # Plot Channel vs Energy
    plt.figure(figsize=(10, 6))
    plt.plot(df['Energy'], df['Channel'], 'b-')
    plt.xlabel('Energy Level (keV)')
    plt.ylabel('Channel')
    plt.title('Channel vs Energy Level')

    # Save the plot as a jpg file
    plot_filename = file_path.split('/')[-1].split('.')[0] + '.jpg'
    plt.savefig(plot_filename, format='jpg')

    # Save data to an Excel file
    excel_filename = file_path.split('/')[-1].split('.')[0] + '.xlsx'
    df.to_excel(excel_filename, index=False)

    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python process_cnf.py <file_path>")
        sys.exit(1)

    process_cnf(sys.argv[1])
