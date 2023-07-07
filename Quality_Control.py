
list=get_subdirectories(s)

for D in Dlist:
    BaseDirList=get_subdirectories(D)

    OdorList = ['ArtCov1/', 'Healthy1k/', 'Healthy100k/', '1propanol/', '2methyl2pentanal/', 'aceticacid/', 'acetone/', 'nonanal/',
                'octanal/', 'pentadecane/', 'ylangylang/']
    ThreshList=[.5,1,10]
    all_dfs = []
    for directory in BaseDirList:

        print(directory)
        print('building dataframe...')
        dfs = [EAG_df_build(os.path.join(directory, odor)) for odor in OdorList]
        all_df = pd.concat(dfs)
        print(all_df.shape)
        all_dfs.append(all_df)

    for directory, all_df in zip(BaseDirList, all_dfs):
        for t in ThreshList:
            print(f"begining quality control for {directory} at threshold of {t}")
            final = FFT_LSTSQ_QC(all_df, t)
            final_T = pd.DataFrame(final.T.dropna(axis=0))

            filename = f"_QC_T_{str(t)}.csv"
            outdir=f"{directory}"
            outdir = outdir.replace('Processed_Data/Official', 'Quality_Controlled_Data')
            print(outdir)
            print(os.path.exists(outdir))

            if os.path.exists(outdir)==False:
                print('making directory...')
                os.makedirs(outdir)
                print('directory made')

            filename = os.path.join(outdir, filename)
            final_T.to_csv(filename)