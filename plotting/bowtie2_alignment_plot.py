import pandas as pd
import altair as alt
from utils import parse_bowtielog


def plot_alignment(bowtie_log: str) -> None:
    """
    Generates alignment plot for the bowtie2 log file
    """

    # read in the numbers
    total_reads, percent_aligned = parse_bowtielog.parse_alignments(bowtie_log)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned

    # Create dataframe
    bowtie_df = pd.DataFrame(
        {"amount": [number_unaligned, number_aligned], "type": ["unaligned", "aligned"]}
    )

    # Generate plot
    plot = (
        alt.Chart(bowtie_df)
        .mark_bar()
        .encode(
            x=alt.X(
                "sum(amount)",
                stack="normalize",
                axis=alt.Axis(format="%"),
                title="Percent",
            ),
            color=alt.Color("type:N", scale=alt.Scale(scheme='dark2')),
            tooltip=[
                alt.Tooltip("amount:Q", title="Number of reads aligned"),
                alt.Tooltip("type:N"),
            ],
        )
        .properties(title="Reads aligned to the human reference genome")
    )
        

    return plot
