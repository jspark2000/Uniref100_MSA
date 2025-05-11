import os
import argparse
import logging
import numpy as np
from evcouplings.align import protocol
from evcouplings.align.alignment import Alignment, read_fasta

logging.basicConfig(
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO
)
logger = logging.getLogger(__name__)


def fetch_sequence(sequence_id, outdir):
    """
    Fetch a protein sequence using UniProt ID
    
    Parameters
    ----------
    sequence_id : str
        UniProt ID of the protein
    outdir : str
        Directory to save the sequence file
        
    Returns
    -------
    str
        Path to the saved sequence file
    """
    os.makedirs(outdir, exist_ok=True)
    sequence_file = os.path.join(outdir, f"{sequence_id}.fasta")
    
    # Use UniProt API to fetch the sequence
    uniprot_url = f'https://www.uniprot.org/uniprot/{sequence_id}.fasta'
    
    if not os.path.exists(sequence_file):
        protocol.fetch_sequence(
            sequence_id, 
            None,  # No input sequence file, will download
            uniprot_url, 
            sequence_file
        )
        logger.info(f"Fetched sequence for {sequence_id}")
    else:
        logger.info(f"Using existing sequence file for {sequence_id}")
    
    return sequence_file


def run_jackhmmer_pipeline(sequence_id, sequence_file, outdir, bitscore=0.3):
    """
    Run the Jackhmmer search pipeline with specified parameters
    
    Parameters
    ----------
    sequence_id : str
        UniProt ID of the protein
    sequence_file : str
        Path to the sequence file
    outdir : str
        Directory to save output files
    bitscore : float
        Bit score threshold in bits per residue
        
    Returns
    -------
    dict
        Output configuration from the pipeline
    """
    # Create output directory
    os.makedirs(outdir, exist_ok=True)
    prefix = os.path.join(outdir, sequence_id)
    
    # Read the sequence to get its length
    with open(sequence_file) as f:
        seq_record = next(read_fasta(f))
    sequence = seq_record[1]
    seq_len = len(sequence)
    
    # Calculate bit score thresholds based on sequence length
    domain_threshold = bitscore * seq_len
    sequence_threshold = domain_threshold
    
    # Set up parameters for jackhmmer search
    kwargs = {
        "prefix": prefix,
        "sequence_id": sequence_id,
        "sequence_file": sequence_file,
        "sequence_download_url": f'https://www.uniprot.org/uniprot/{sequence_id}.fasta',
        "region": None,
        "first_index": 1,
        "use_bitscores": True,
        "domain_threshold": domain_threshold,
        "sequence_threshold": sequence_threshold,
        "database": "database",
        "database_path": "/path/to/uniref100.fasta",
        "iterations": 5,
        "cpu": 32,
        "nobias": False,
        "reuse_alignment": False,
        "checkpoints_hmm": False,
        "checkpoints_ali": False,
        "jackhmmer": "/data01/genbiolab/modules/anaconda3/2024.10/envs/evecp/bin/jackhmmer",
        "extract_annotation": True
    }
    
    kwargs["database"] = "database_path"
    
    logger.info(f"Running jackhmmer search for {sequence_id} with bitscore {bitscore}")
    
    search_cfg = protocol.jackhmmer_search(**kwargs)
    
    return search_cfg


def filter_alignment(alignment_file, target_seq_id, outdir, 
                     min_coverage=0.7, column_occupancy=0.7):
    """
    Filter the alignment based on coverage and column occupancy criteria
    
    Parameters
    ----------
    alignment_file : str
        Path to the raw alignment file
    target_seq_id : str
        ID of the target sequence
    outdir : str
        Directory to save filtered alignment
    min_coverage : float
        Minimum coverage of target sequence (default: 0.7)
    column_occupancy : float
        Minimum column occupancy (default: 0.7)
        
    Returns
    -------
    tuple
        (filtered_alignment_file, stats_dict) with path to filtered alignment 
        and alignment statistics
    """
    # Read the alignment
    with open(alignment_file) as f:
        ali = Alignment.from_file(f, format="stockholm")
    
    # Find the index of the target sequence
    target_seq_index = 0
    if target_seq_id not in ali.ids[0]:
        for i, seq_id in enumerate(ali.ids):
            if target_seq_id in seq_id:
                target_seq_index = i
                break
    
    # Calculate coverage for each sequence
    target_seq = ali.matrix[target_seq_index]
    seq_len = (target_seq != "-").sum()
    
    # Filter sequences that align to at least 50% of target protein
    coverage_mask = np.zeros(len(ali.ids), dtype=bool)
    for i, seq in enumerate(ali.matrix):
        # Calculate coverage (non-gap positions aligned to target sequence)
        aligned_positions = ((seq != "-") & (target_seq != "-")).sum()
        coverage = aligned_positions / seq_len
        coverage_mask[i] = coverage >= 0.5
    
    # Keep the target sequence regardless of coverage
    coverage_mask[target_seq_index] = True
    
    # Apply sequence filter
    filtered_ids = [ali.ids[i] for i in range(len(ali.ids)) if coverage_mask[i]]
    filtered_matrix = ali.matrix[coverage_mask]
    
    # Filter columns with at least 70% residue occupancy
    column_mask = np.zeros(ali.matrix.shape[1], dtype=bool)
    for j in range(ali.matrix.shape[1]):
        occupancy = (filtered_matrix[:, j] != "-").sum() / filtered_matrix.shape[0]
        column_mask[j] = occupancy >= column_occupancy
    
    # Apply column filter, but ensure we keep columns where target has residues
    target_residue_mask = target_seq != "-"
    final_column_mask = column_mask | target_residue_mask
    
    # Create filtered alignment
    final_matrix = filtered_matrix[:, final_column_mask]
    filtered_ali = Alignment(final_matrix, filtered_ids)
    
    # Save filtered alignment
    filtered_alignment_file = os.path.join(outdir, f"{target_seq_id}_filtered.a2m")
    with open(filtered_alignment_file, "w") as f:
        filtered_ali.write(f, format="fasta")
    
    # Calculate statistics
    n_seqs = len(filtered_ids)
    coverage_columns = final_column_mask.sum() / len(target_residue_mask)
    
    stats = {
        "num_sequences": n_seqs,
        "target_seq_length": seq_len,
        "coverage": coverage_columns,
        "filtered_columns": final_column_mask.sum(),
        "original_columns": len(target_residue_mask)
    }
    
    logger.info(f"Filtered alignment: {n_seqs} sequences, coverage: {coverage_columns:.2f}")
    
    return filtered_alignment_file, stats


def evaluate_alignment_quality(stats, target_seq_len):
    """
    Evaluate alignment quality based on coverage and sequence count criteria
    
    Parameters
    ----------
    stats : dict
        Alignment statistics
    target_seq_len : int
        Length of the target sequence
        
    Returns
    -------
    tuple
        (quality_score, meets_criteria) with a numeric score for quality
        and boolean indicating if it meets all criteria
    """
    coverage = stats["coverage"]
    num_sequences = stats["num_sequences"]
    
    # Calculate coverage relative to target sequence length
    l_cov = coverage * target_seq_len
    
    # Check primary criteria: Lcov ≥ 0.8L and 100,000 ≥ N ≥ 10L
    primary_coverage = l_cov >= 0.8 * target_seq_len
    primary_seq_count = 10 * target_seq_len <= num_sequences <= 100000
    
    # Check secondary criteria: Lcov ≥ 0.7L and N ≤ 200,000
    secondary_coverage = l_cov >= 0.7 * target_seq_len
    secondary_seq_count = num_sequences <= 200000
    
    # Calculate a quality score (higher is better)
    # Favor alignments with higher coverage and appropriate sequence count
    quality_score = coverage * min(1.0, num_sequences / (10 * target_seq_len))
    
    # Penalize having too many sequences
    if num_sequences > 100000:
        quality_score *= (100000 / num_sequences)
    
    # Check if criteria are met
    meets_primary = primary_coverage and primary_seq_count
    meets_secondary = secondary_coverage and secondary_seq_count
    meets_criteria = meets_primary or meets_secondary or num_sequences > 0
    
    return quality_score, meets_criteria


def select_best_alignment(sequence_id, sequence_file, outdir, 
                          bitscores=None):
    """
    Run jackhmmer with different bit scores and select the best alignment
    
    Parameters
    ----------
    sequence_id : str
        UniProt ID of the protein
    sequence_file : str
        Path to the sequence file
    outdir : str
        Directory to save output files
    bitscores : list, optional
        List of bit score thresholds to try
        
    Returns
    -------
    tuple
        (best_alignment_file, best_stats) with path to best alignment and its statistics
    """
    if bitscores is None:
        bitscores = [0.2, 0.25, 0.3, 0.35, 0.4, 0.5]
    
    # Read target sequence to get its length
    with open(sequence_file) as f:
        seq_record = next(read_fasta(f))
    target_seq_len = len(seq_record[1])
    
    best_score = -1
    best_alignment = None
    best_stats = None
    
    # Try each bit score threshold
    for bitscore in bitscores:
        bitscore_dir = os.path.join(outdir, f"bitscore_{bitscore}")
        os.makedirs(bitscore_dir, exist_ok=True)
        
        # Run jackhmmer search
        search_cfg = run_jackhmmer_pipeline(
            sequence_id, sequence_file, bitscore_dir, bitscore
        )
        
        # Filter alignment
        filtered_alignment, stats = filter_alignment(
            search_cfg["raw_alignment_file"], sequence_id, bitscore_dir
        )
        
        # Evaluate alignment quality
        quality_score, meets_criteria = evaluate_alignment_quality(stats, target_seq_len)
        
        logger.info(f"Bitscore {bitscore}: quality score = {quality_score:.4f}, "
                   f"meets criteria = {meets_criteria}")
        
        # Update best alignment if this one is better
        if meets_criteria and quality_score > best_score:
            best_score = quality_score
            best_alignment = filtered_alignment
            best_stats = stats
    
    if best_alignment is None:
        logger.warning("No alignment meets criteria, selecting based on sequence count")
        best_score = -1
        for bitscore in bitscores:
            bitscore_dir = os.path.join(outdir, f"bitscore_{bitscore}")
            filtered_alignment = os.path.join(bitscore_dir, f"{sequence_id}_filtered.a2m")
            
            if os.path.exists(filtered_alignment):
                with open(filtered_alignment) as f:
                    ali = Alignment.from_file(f, format="fasta")
                
                if len(ali.ids) > best_score:
                    best_score = len(ali.ids)
                    best_alignment = filtered_alignment
                    
                    stats = {
                        "num_sequences": len(ali.ids),
                        "coverage": 0,
                        "filtered_columns": ali.L,
                        "original_columns": ali.L
                    }
                    best_stats = stats
    
    # Copy the best alignment to the main output directory
    if best_alignment:
        final_alignment = os.path.join(outdir, f"{sequence_id}_final.a2m")
        with open(best_alignment) as f_in:
            with open(final_alignment, "w") as f_out:
                f_out.write(f_in.read())
        
        logger.info(f"Selected best alignment with {best_stats['num_sequences']} sequences "
                   f"and coverage {best_stats.get('coverage', 0):.2f}")
        
        return final_alignment, best_stats
    else:
        logger.error("Failed to generate any valid alignments")
        return None, None


def main():
    parser = argparse.ArgumentParser(description="Generate MSAs for protein families")
    
    parser.add_argument("--sequence_id", required=True, help="UniProt ID of the protein")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--database", default="/path/to/uniref100.fasta", 
                       help="Path to UniRef100 database")
    parser.add_argument("--cpu", type=int, default=4, 
                       help="Number of CPUs to use")
    parser.add_argument("--bitscores", type=float, nargs="+", 
                       default=[0.2, 0.25, 0.3, 0.35, 0.4, 0.5],
                       help="Bit score thresholds to try")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Step 1: Fetch the protein sequence
    sequence_file = fetch_sequence(args.sequence_id, args.outdir)
    
    # Step 2: Run the pipeline with different bit scores and select the best alignment
    best_alignment, stats = select_best_alignment(
        args.sequence_id, sequence_file, args.outdir, args.bitscores
    )
    
    if best_alignment:
        # Save alignment statistics
        stats_file = os.path.join(args.outdir, f"{args.sequence_id}_stats.txt")
        with open(stats_file, "w") as f:
            for key, value in stats.items():
                f.write(f"{key}: {value}\n")
        
        logger.info(f"MSA generation complete. Output saved to {best_alignment}")
    else:
        logger.error("Failed to generate MSA")


if __name__ == "__main__":
    main() 