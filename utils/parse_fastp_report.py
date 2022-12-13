from bs4 import BeautifulSoup
import pandas as pd

def parse_fastp(fastp_report: str) -> pd.DataFrame:
    """
    Parses the relevant information about reads from fastp report
    """
    
    summary_information = {}

    with open(fastp_report, "r") as f:
        soup = BeautifulSoup(f, "html.parser")

    for table in soup.find_all("table", class_="summary_table")[:4]:
        for row in table.find_all("tr"):
            cells = row.find_all("td")
            key = cells[0].text
            value = cells[1].text
            summary_information[key] = value
            
    df = pd.DataFrame.from_dict(summary_information, orient="index", columns=["value"])
    df = df.rename_axis("description").reset_index()
    return df