import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MultipleLocator
from numpy.polynomial.polynomial import polyfit
from pathlib import Path

def analyze_em27_dataset(
        filepath,
        filepath_corr_XCO2,
        filepath_corr_XCH4,
        station,
        date,
        filepath_corr_XAIR=None,
        save_dir=None,
        ylim_XAIR_before=0.01,
        ylim_XCO2_before=5,
        ylim_XCH4_before=0.03,
        ylim_XAIR_after=0.01, 
        ylim_XCO2_after=2,
        ylim_XCH4_after=0.01,
        start = None ,
        end = None , 
        fit_summary_file = None , 
        Random_values = False,
        save_XAIR = False
        ):
    
    ########################
    #Plots before correction
    ########################
    
    #Read file before the correction
    ds = pd.read_csv(filepath,skipinitialspace=True)

    #Incase the data has random values 
    if Random_values :
            ds =ds [(ds['XCO2'] > 400) & 
                (ds['XCO2'] <450) &
                (ds['XCH4'] > 1.5 ) &
                (ds['XCH4'] < 2) &
                (ds['XAIR'] > 0.95 ) &
                (ds['XAIR'] < 1.05 )
                ]

    #Get the values we need from ds
    ds_small = ds[['UTC','LocalTime','JulianDate','appSZA','XAIR','XCO2','XCH4']]

    #Filter out invalid values
    ds_small_filt = ds_small[ds_small['appSZA']<= 86]

    #compute the airmass - before the correction
    airmass = 1/np.cos(np.deg2rad(ds_small_filt['appSZA']))

    #Get the real dates from julian date
    time = pd.to_datetime(ds_small_filt['JulianDate'],unit='D' , origin='julian')

    # Add time into dataframe so we can filter by date
    ds_small_filt = ds_small_filt.copy()
    ds_small_filt['time'] = time

    # Optional date filtering
    if start is not None:
        ds_small_filt = ds_small_filt[ds_small_filt['time'] >= pd.to_datetime(start)]

    if end is not None:
        ds_small_filt = ds_small_filt[ds_small_filt['time'] <= pd.to_datetime(end)]

    time = ds_small_filt['time']

    #compute the mean values before the correction
    mean_XCO2_before = np.mean(ds_small_filt['XCO2'])
    mean_XCH4_before = np.mean(ds_small_filt['XCH4'])
    mean_XAIR_before = np.mean(ds_small_filt['XAIR'])

    #Plot XCO2 vs Date
    fig,ax = plt.subplots(figsize=(10,5))
    ax.scatter(time , ds_small_filt['XCO2'], s=2 ,color='blue', label='XCO2')
    ax.set_xlabel('Date')
    ax.set_ylim(mean_XCO2_before - ylim_XCO2_before, mean_XCO2_before + ylim_XCO2_before)
    ax.set_ylabel('XCO2 (ppm)')
    ax.set_title(f'XCO2 over Time - {station}-({date}) Before Correction')
    if save_dir is not None:
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XCO2_before-Time_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    plt.show()

    #Plot XCO2 vs Airmass/appSZA before correction

    fig,ax = plt.subplots(figsize=(15,4))
    x_XCO2_before = ds_small_filt['appSZA']
    y_XCO2_before = ds_small_filt['XCO2']
    ax.scatter(x_XCO2_before,y_XCO2_before,s=2 , color='blue')

    def sza_to_airmass(sza_deg):
        return 1/np.cos(np.deg2rad(sza_deg))

    def airmass_to_sza(airmass):
        airmass = np.asarray(airmass)
        airmass = np.clip(airmass, 1, None)
        return np.rad2deg(np.arccos(1/airmass))
    
    ax.margins(x=0)
    ax.set_xlim(x_XCO2_before.min(),86)
    
    secax = ax.secondary_xaxis('top', functions=(sza_to_airmass,airmass_to_sza))
    secax.set_xlabel('Airmass')
    ax.set_ylim(mean_XCO2_before - ylim_XCO2_before, mean_XCO2_before + ylim_XCO2_before)
    ax.axhline(y=np.nanmean(y_XCO2_before), color='black',label='Mean XCO2(ppm)')
    secax.set_xticks(np.arange(1,10,1))
    secax.xaxis.set_minor_locator(MultipleLocator(0.5))
    secax.tick_params(axis='x', which='minor', length=3)
    secax.tick_params(axis='x',which='major',length=6)
    ax.set_xlabel('appSZA_deg')
    ax.set_ylabel('XCO2 (ppm)')
    ax.set_title(f'XCO2 vs appSZA/airmass - {station}-({date})-Before Correction , Mean +- {ylim_XCO2_before} ppm')
    ax.legend()
    if save_dir is not None:
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XCO2_before-SZA-Airmass_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    
    plt.show()

    #Plot XCH4 vs Date
    fig,ax = plt.subplots(figsize=(10,5))
    ax.scatter(time , ds_small_filt['XCH4'], s=2 ,color='green', label='XCH4')
    ax.set_xlabel('Date')
    ax.set_ylim(mean_XCH4_before - ylim_XCH4_before, mean_XCH4_before + ylim_XCH4_before)
    ax.set_ylabel('XCH4 (ppm)')
    ax.set_title(f'XCH4 over Time - {station}-({date}) Before Correction')
    if save_dir is not None:
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XCH4_before-Time_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    
    plt.show()

    #Plot XCH4 vs Airmass/appSZA before correction

    fig,ax = plt.subplots(figsize=(15,4))
    x_XCH4_before= ds_small_filt['appSZA']
    y_XCH4_before= ds_small_filt['XCH4']
    ax.scatter(x_XCH4_before,y_XCH4_before,s=2,color='green')
    ax.margins(x=0)
    ax.set_xlim(x_XCH4_before.min(),86)
    secax = ax.secondary_xaxis('top', functions=(sza_to_airmass,airmass_to_sza))
    secax.set_xlabel('Airmass')
    ax.set_ylim(mean_XCH4_before - ylim_XCH4_before, mean_XCH4_before + ylim_XCH4_before)
    secax.set_xticks(np.arange(1,10,1))
    secax.xaxis.set_minor_locator(MultipleLocator(0.5))
    secax.tick_params(axis='x', which='minor', length=3)
    secax.tick_params(axis='x',which='major',length=6)
    ax.axhline(y=np.nanmean(y_XCH4_before), color='black',label='Mean XCH4(ppm)')
    ax.set_xlabel('appSZA_deg')
    ax.set_ylabel('XCH4 (ppm)')
    ax.set_title(f'XCH4 vs appSZA/airmass - {station}-({date})-Before Correction , Mean +- {ylim_XCH4_before} ppm')
    ax.legend()
    if save_dir is not None:
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XCH4_before-SZA-Airmass_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    
    
    plt.show()

    #Plot XAIR vs Date
    fig,ax = plt.subplots(figsize=(10,5))
    ax.scatter(time , ds_small_filt['XAIR'], s=2 ,color='red', label='XAIR')
    ax.set_xlabel('Date')
    ax.set_ylim(mean_XAIR_before - ylim_XAIR_before, mean_XAIR_before + ylim_XAIR_before)
    ax.set_ylabel('XAIR ')
    ax.set_title(f'XAiR over Time - {station}-({date}) Before Correction')
    if save_dir is not None and save_XAIR :
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XAIR_before-Time_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    
    plt.show()

    #Plot XAIR vs Airmass/appSZA before correction

    fig,ax = plt.subplots(figsize=(15,4))
    x_XAIR_before= ds_small_filt['appSZA']
    y_XAIR_before= ds_small_filt['XAIR']
    ax.scatter(x_XAIR_before,y_XAIR_before,s=2,color='red')
    ax.margins(x=0)
    ax.set_xlim(x_XAIR_before.min(),86)
    secax = ax.secondary_xaxis('top', functions=(sza_to_airmass,airmass_to_sza))
    secax.set_xlabel('Airmass')
    ax.set_ylim(mean_XAIR_before - ylim_XAIR_before, mean_XAIR_before + ylim_XAIR_before)
    secax.set_xticks(np.arange(1,10,1))
    secax.xaxis.set_minor_locator(MultipleLocator(0.5))
    secax.tick_params(axis='x', which='minor', length=3)
    secax.tick_params(axis='x',which='major',length=6)
    ax.axhline(y=np.nanmean(y_XAIR_before), color='black',label='Mean XAIR)')
    ax.set_xlabel('appSZA_deg')
    ax.set_ylabel('XAIR ')
    ax.set_title(f'XAIR vs appSZA/airmass - {station}-({date})-Before Correction , Mean +- {ylim_XAIR_before}')
    ax.legend()
    if save_dir is not None and save_XAIR :
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XAIR_before-SZA-Airmass_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    
    
    plt.show()

    #######################
    #Plots after correction
    #######################

    #Get the corrected dataset
    ds_corr_XCO2 = pd.read_csv(filepath_corr_XCO2,names = ['JulianDate','appSZA','XCO2'])

    #Incase the data has random values 
    if Random_values :
            ds_corr_XCO2 = ds_corr_XCO2 [
                (ds_corr_XCO2['XCO2'] >400) & 
                (ds_corr_XCO2['XCO2'] <450) 
                ]

    #Filter the data for appsza<86
    ds_corr_XCO2_filt = ds_corr_XCO2[ds_corr_XCO2['appSZA']<= 86]

    #Get time for corrected XCO2
    time_XCO2 = pd.to_datetime(ds_corr_XCO2_filt['JulianDate'],unit='D' , origin='julian')

    ds_corr_XCO2_filt = ds_corr_XCO2_filt.copy()
    ds_corr_XCO2_filt['time'] = time_XCO2

    #Optional date filtering
    if start is not None:
        ds_corr_XCO2_filt = ds_corr_XCO2_filt[ds_corr_XCO2_filt['time'] >= pd.to_datetime(start)]

    if end is not None:
        ds_corr_XCO2_filt = ds_corr_XCO2_filt[ds_corr_XCO2_filt['time'] <= pd.to_datetime(end)]

    time_XCO2 = ds_corr_XCO2_filt['time']

    #find airmass for corrected XCO2 
    airmass_corr = 1/np.cos(np.deg2rad(ds_corr_XCO2_filt['appSZA']))

    #compute mean column concentration for corrected XCO2
    mean_XCO2_corr = np.mean(ds_corr_XCO2_filt['XCO2'])

    #Plot XCO2 vs Date
    fig,ax = plt.subplots(figsize=(10,5))
    ax.scatter(time_XCO2 , ds_corr_XCO2_filt['XCO2'], s=2 ,color='blue', label='XCO2')
    ax.set_xlabel('Date')
    ax.set_ylim(mean_XCO2_corr - ylim_XCO2_after, mean_XCO2_corr + ylim_XCO2_after)
    ax.set_ylabel('XCO2 (ppm)')
    ax.set_title(f'XCO2 over Time - {station}-({date}) After Correction')
    if save_dir is not None:
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XCO2_after-Time_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    
    
    plt.show()
    
    #Plot corrected XCO2 vs Airmass/appSZA 
    fig,ax = plt.subplots(figsize=(15,4))
    x_XCO2_after = ds_corr_XCO2_filt['appSZA']
    y_XCO2_after = ds_corr_XCO2_filt['XCO2']
    ax.scatter(x_XCO2_after,y_XCO2_after,s=2 , color='blue')

    #Add regression lines for corrected XCO2 with SZA<=70
    fit_mask = x_XCO2_after <=70

    x_fit_data = x_XCO2_after[fit_mask]
    y_fit_data = y_XCO2_after[fit_mask]

    x_fit = np.linspace(x_fit_data.min(), 70, 500)

    # 1st degree regression
    c1, b1 = polyfit(x_fit_data, y_fit_data, 1)
    y_fit_1 = c1 + b1*x_fit

    ax.plot(
    x_fit,
    y_fit_1,
    color='purple',
    label=f"Linear fit ≤70°: {b1:.2e}x + {c1:.2f}"
    )

    # 2nd degree regression
    c2, b2, a2 = polyfit(x_fit_data, y_fit_data, 2)
    y_fit_2 = c2 + b2*x_fit + a2*x_fit**2

    ax.plot(
    x_fit,
    y_fit_2,
    color='cyan',
    label=f"Quadratic fit ≤70°: {a2:.2e}x$^2$ + {b2:.2e}x + {c2:.2f}"
    )

    # Save fit results automatically
    if fit_summary_file is not None:

        fit_summary_file = Path(fit_summary_file)
        
        gas_name = "XCO2"

        fit_row = pd.DataFrame([{
            "station": station,
            "date": date,
            "gas": gas_name ,
            "linear_slope": round(b1, 8),
            "quadratic_curvature": round(a2, 10),
            "mean": round(mean_XCO2_corr, 3),
            "n_points": len(x_fit_data)
        }])

        if fit_summary_file.exists():
            old = pd.read_csv(fit_summary_file, sep=";", decimal=",")
            # Remove old row with same station/date/gas
            old = old[
                ~(
                    (old["station"] == station) &
                    (old["date"] == date) &
                    (old["gas"] == gas_name)
                )
        ]
            out = pd.concat([old, fit_row], ignore_index=True)
        else:
            out = fit_row

        # sort alphabetically by station, then by date, then by gas
        out = out.sort_values(
        by=["station", "date", "gas"],
        ascending=[True, True, True]
        )

        out.to_csv(fit_summary_file, sep=";", decimal=",", index=False)

    ax.margins(x=0)
    ax.set_xlim(x_XCO2_after.min(),86)

    secax = ax.secondary_xaxis('top', functions=(sza_to_airmass,airmass_to_sza))
    secax.set_xlabel('Airmass')
    ax.set_ylim(mean_XCO2_corr - ylim_XCO2_after, mean_XCO2_corr + ylim_XCO2_after)
    secax.set_xticks(np.arange(1,10,1))
    secax.xaxis.set_minor_locator(MultipleLocator(0.5))
    secax.tick_params(axis='x', which='minor', length=3)
    secax.tick_params(axis='x',which='major',length=6)
    ax.axhline(y=np.nanmean(y_XCO2_after), color='black',label='Mean XCO2(ppm)')
    ax.set_xlabel('appSZA_deg')
    ax.set_ylabel('XCO2 (ppm)')
    ax.set_title(f'XCO2 vs appSZA/airmass - {station}-({date}) After Correction , Mean +- {ylim_XCO2_after} ppm')
    plt.legend()
    if save_dir is not None:
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XCO2_after-SZA-Airmass_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')

    
    plt.show()

    #Get corrected XCH4 dataset
    ds_corr_XCH4 = pd.read_csv(filepath_corr_XCH4 ,names = ['JulianDate','appSZA','XCH4'])

    #Incase the data has random values 
    if Random_values :
            ds_corr_XCH4 = ds_corr_XCH4 [
                (ds_corr_XCH4['XCH4'] > 1.7) & 
                (ds_corr_XCH4['XCH4'] < 2) 
                ]

    #filter XCH4  data
    ds_corr_XCH4_filt = ds_corr_XCH4[ds_corr_XCH4['appSZA'] <= 86 ]

    #Get time for XCH4 
    time_XCH4 = pd.to_datetime(ds_corr_XCH4_filt['JulianDate'] , unit= 'D' , origin='julian')

    ds_corr_XCH4_filt = ds_corr_XCH4_filt.copy()
    ds_corr_XCH4_filt['time'] = time_XCH4

    #Optional date filtering 
    if start is not None:
        ds_corr_XCH4_filt = ds_corr_XCH4_filt[ds_corr_XCH4_filt['time'] >= pd.to_datetime(start)]

    if end is not None:
        ds_corr_XCH4_filt = ds_corr_XCH4_filt[ds_corr_XCH4_filt['time'] <= pd.to_datetime(end)]

    time_XCH4 = ds_corr_XCH4_filt['time']

    #airmass for corrected XCH4 
    airmass_corr_XCH4 = 1/np.cos(np.deg2rad(ds_corr_XCH4_filt['appSZA']))

    #mean 
    mean_XCH4_after = np.mean(ds_corr_XCH4_filt['XCH4'])

    #Plot corrected XCH4 vs Date
    fig,ax = plt.subplots(figsize=(10,5))
    ax.scatter(time_XCH4 , ds_corr_XCH4_filt['XCH4'], s=2 ,color='green', label='XCH4')
    ax.set_xlabel('Date')
    ax.set_ylim(mean_XCH4_after - ylim_XCH4_after, mean_XCH4_after + ylim_XCH4_after)
    ax.set_ylabel('XCH4 (ppm)')
    ax.set_title(f'XCH4 over Time - {station}-({date}) After Correction')
    if save_dir is not None:
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XCH4_after-Time_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    plt.show()

    #Plot XCH4 vs Airmass/appSZA after correction

    fig,ax = plt.subplots(figsize=(15,4))
    x_XCH4_after= ds_corr_XCH4_filt['appSZA']
    y_XCH4_after= ds_corr_XCH4_filt['XCH4']
    ax.scatter(x_XCH4_after,y_XCH4_after,s=2,color='green')

    #Add regression lines for corrected XCH4 with SZA<=70
    fit_mask = x_XCH4_after <=70

    x_fit_data = x_XCH4_after[fit_mask]
    y_fit_data = y_XCH4_after[fit_mask]

    x_fit = np.linspace(x_fit_data.min(), 70, 500)

    # 1st degree regression
    c1, b1 = polyfit(x_fit_data, y_fit_data, 1)
    y_fit_1 = c1 + b1*x_fit

    ax.plot(
    x_fit,
    y_fit_1,
    color='purple',
    label=f"Linear fit ≤70°: {b1:.2e}x + {c1:.2f}"
    )

    # 2nd degree regression
    c2, b2, a2 = polyfit(x_fit_data, y_fit_data, 2)
    y_fit_2 = c2 + b2*x_fit + a2*x_fit**2

    ax.plot(
    x_fit,
    y_fit_2,
    color='cyan',
    label=f"Quadratic fit ≤70°: {a2:.2e}x$^2$ + {b2:.2e}x + {c2:.2f}"
    )

    # Save fit results automatically
    if fit_summary_file is not None:

        fit_summary_file = Path(fit_summary_file)

        gas_name = "XCH4"
        
        fit_row = pd.DataFrame([{
            "station": station,
            "date": date,
            "gas":gas_name,
            "linear_slope": round(b1, 8),
            "quadratic_curvature": round(a2, 10),
            "mean": round(mean_XCH4_after, 3),
            "n_points": len(x_fit_data)
        }])

        if fit_summary_file.exists():
            old = pd.read_csv(fit_summary_file, sep=";", decimal=",")
            # Remove old row with same station/date/gas
            old = old[
                ~(
                    (old["station"] == station) &
                    (old["date"] == date) &
                    (old["gas"] == gas_name)
                )
        ]
            out = pd.concat([old, fit_row], ignore_index=True)
        else:
            out = fit_row

        # sort alphabetically by station, then by date, then by gas
        out = out.sort_values(
        by=["station", "date", "gas"],
        ascending=[True, True, True]
        )

        out.to_csv(fit_summary_file, sep=";", decimal=",", index=False)

    ax.margins(x=0)
    ax.set_xlim(x_XCH4_after.min(),86)
    secax = ax.secondary_xaxis('top', functions=(sza_to_airmass,airmass_to_sza))
    secax.set_xlabel('Airmass')
    ax.set_ylim(mean_XCH4_after - ylim_XCH4_after, mean_XCH4_after + ylim_XCH4_after)
    secax.set_xticks(np.arange(1,10,1))
    secax.xaxis.set_minor_locator(MultipleLocator(0.5))
    secax.tick_params(axis='x', which='minor', length=3)
    secax.tick_params(axis='x',which='major',length=6)
    ax.axhline(y=np.nanmean(y_XCH4_after), color='black',label='Mean XCH4(ppm)')
    ax.set_xlabel('appSZA_deg')
    ax.set_ylabel('XCH4 (ppm)')
    ax.set_title(f'XCH4 vs appSZA/airmass - {station}-({date})-After Correction , Mean +- {ylim_XCH4_after} ppm')
    ax.legend()
    if save_dir is not None:
        output_path = Path(save_dir) / station / date
        output_path.mkdir(parents=True,exist_ok=True)
        fig.savefig(output_path / f"XCH4_after-SZA-Airmass_{station}-({date}).png",
                    dpi=100,
                    bbox_inches='tight')
    
    plt.show()

    if filepath_corr_XAIR is not None: 

        #XAIR corrected dataset
        ds_corr_XAIR = pd.read_csv(filepath_corr_XAIR ,names = ['JulianDate','appSZA','XAIR'])

        #Incase the data has random values 
        if Random_values :
            ds_corr_XAIR = ds_corr_XAIR [
            (ds_corr_XAIR['XAIR'] > 0.98) & 
            (ds_corr_XAIR['XAIR'] < 1.02) 
            ]

        #filter XAIR  data
        ds_corr_XAIR_filt = ds_corr_XAIR[ds_corr_XAIR['appSZA'] <= 86 ]

        #Get time for XAIR 
        time_XAIR = pd.to_datetime(ds_corr_XAIR_filt['JulianDate'] , unit= 'D' , origin='julian')

        # Add time into dataframe so we can filter by date
        ds_corr_XAIR_filt = ds_corr_XAIR_filt.copy()
        ds_corr_XAIR_filt['time'] = time_XAIR

        #Optional date filtering
        if start is not None:
            ds_corr_XAIR_filt = ds_corr_XAIR_filt[ds_corr_XAIR_filt['time'] >= pd.to_datetime(start)]

        if end is not None:
            ds_corr_XAIR_filt = ds_corr_XAIR_filt[ds_corr_XAIR_filt['time'] <= pd.to_datetime(end)]

        time_XAIR = ds_corr_XAIR_filt['time']

        #airmass for corrected XAIR 
        airmass_corr_XAIR = 1/np.cos(np.deg2rad(ds_corr_XAIR_filt['appSZA']))

        #mean 
        mean_XAIR_after = np.mean(ds_corr_XAIR_filt['XAIR'])

        #Plot XAIR vs Airmass/appSZA after correction

        fig,ax = plt.subplots(figsize=(15,4))
        x_XAIR_after= ds_corr_XAIR_filt['appSZA']
        y_XAIR_after= ds_corr_XAIR_filt['XAIR']
        ax.scatter(x_XAIR_after,y_XAIR_after,s=2,color='red')

        #Add regression lines for corrected XAIR with SZA<=70
        fit_mask = x_XAIR_after <=70

        x_fit_data = x_XAIR_after[fit_mask]
        y_fit_data = y_XAIR_after[fit_mask]

        x_fit = np.linspace(x_fit_data.min(), 70, 500)

        # 1st degree regression
        c1, b1 = polyfit(x_fit_data, y_fit_data, 1)
        y_fit_1 = c1 + b1*x_fit

        ax.plot(
        x_fit,
        y_fit_1,
        color='purple',
        label=f"Linear fit ≤70°: {b1:.2e}x + {c1:.2f}"
        )

        # 2nd degree regression
        c2, b2, a2 = polyfit(x_fit_data, y_fit_data, 2)
        y_fit_2 = c2 + b2*x_fit + a2*x_fit**2

        ax.plot(
        x_fit,
        y_fit_2,
        color='cyan',
        label=f"Quadratic fit ≤70°: {a2:.2e}x$^2$ + {b2:.2e}x + {c2:.2f}")

         # Save fit results automatically
        if fit_summary_file is not None:

            fit_summary_file = Path(fit_summary_file)

            gas_name="XAIR"
            
            fit_row = pd.DataFrame([{
                "station": station,
                "date": date,
                "gas": gas_name,
                "linear_slope": round(b1, 8),
                "quadratic_curvature": round(a2, 10),
                "mean": round(mean_XAIR_after, 3),
                "n_points": len(x_fit_data)
            }])

        if fit_summary_file.exists():
            old = pd.read_csv(fit_summary_file, sep=";", decimal=",")
            # Remove old row with same station/date/gas
            old = old[
                ~(
                    (old["station"] == station) &
                    (old["date"] == date) &
                    (old["gas"] == gas_name)
                )
        ]
            out = pd.concat([old, fit_row], ignore_index=True)
        else:
            out = fit_row

        # sort alphabetically by station, then by date, then by gas
        out = out.sort_values(
        by=["station", "date", "gas"],
        ascending=[True, True, True]
        )

        out.to_csv(fit_summary_file, sep=";", decimal=",", index=False)

        ax.margins(x=0)
        ax.set_xlim(x_XAIR_after.min(),86)
        secax = ax.secondary_xaxis('top', functions=(sza_to_airmass,airmass_to_sza))
        secax.set_xlabel('Airmass')
        ax.set_ylim(mean_XAIR_after - ylim_XAIR_after, mean_XAIR_after + ylim_XAIR_after)
        secax.set_xticks(np.arange(1,10,1))
        secax.xaxis.set_minor_locator(MultipleLocator(0.5))
        secax.tick_params(axis='x', which='minor', length=3)
        secax.tick_params(axis='x',which='major',length=6)
        ax.axhline(y=np.nanmean(y_XAIR_after), color='black',label='Mean XAIR')
        ax.set_xlabel('appSZA_deg')
        ax.set_ylabel('XAIR')
        ax.set_title(f'XAIR vs appSZA/airmass - {station}-({date})-After Correction , Mean +- {ylim_XAIR_after}')
        ax.legend()
        if save_dir is not None:
            output_path = Path(save_dir) / station / date
            output_path.mkdir(parents=True,exist_ok=True)
            fig.savefig(output_path / f"XAIR_after-SZA-Airmass_{station}-({date}).png",
                        dpi=100,
                        bbox_inches='tight')
        
        plt.show()
    