import datetime
import pandas as pd
import sys
import os


def make_periods_csv(first_day, last_day, stepsize_days, save_dir):
    """
    Create a csv file containing start and end dates for Kalman filter update periods.
    For example, setting first_day='20180501', last_day='20180601', stepsize_days=7
    will generate a csv with weekly start/end days: 20180501-20180508, 20180508-20180515, etc.

    Arguments
        first_day     [str] : First day of Kalman filter inversion period (yyymmdd)
        last_day      [str] : Last day of Kalman filter inversion period (yyymmdd)
        stepsize_days [int] : Update frequency (number of days)
        save_dir      [str] : Where to save the periods.csv file
    """

    # Initialize empty lists for start and end days
    starts = []
    ends = []

    dt = datetime.datetime.strptime(first_day, "%Y%m%d")
    dt_max = datetime.datetime.strptime(last_day, "%Y%m%d")
    delta = datetime.timedelta(days=stepsize_days)
    remainder = (dt - dt_max).days % stepsize_days

    # Validate entire period is divisible by step size
    if dt + delta > dt_max:
        raise Exception(
            "Stepsize exceeds the end date for first period. Please reduce UpdateFreqDays."
        )
    elif remainder != 0:
        print(
            "Warning: Number of days between start date and end"
            + f" date is not divisible by UpdateFreqDays: {stepsize_days}."
            + f" \nConsider updating end date to have +{remainder} days"
            + f" or reducing by -{stepsize_days - remainder} days."
        )

    # Use while loop to populate starts and ends
    while dt < dt_max:
        dt_start_str = str(dt)[0:10].replace("-", "")
        dt_start_int = int(dt_start_str)
        delta = datetime.timedelta(days=stepsize_days)
        dt += delta
        dt_end_str = str(dt)[0:10].replace("-", "")
        dt_end_int = int(dt_end_str)
        if dt <= dt_max:
            starts.append(dt_start_int)
            ends.append(dt_end_int)

    df = pd.DataFrame(list(zip(starts, ends)), columns=["Starts", "Ends"])
    df["period_number"] = df.index.values + 1
    df.to_csv(os.path.join(save_dir, "periods.csv"), header=True, index=False)


if __name__ == "__main__":
    # Inputs
    first_day = sys.argv[1]
    last_day = sys.argv[2]
    stepsize_days = sys.argv[3]
    save_dir = sys.argv[4]

    # Run the script
    make_periods_csv(first_day, last_day, int(stepsize_days), save_dir)
