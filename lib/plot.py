def pca_exprs(df, ntop:int = 5000, log2 = True, pseudocount = 0.001, n_components:int = 2, metadata = None):
    if log2:
        df = np.log2(df + pseudocount)
    top = df.var(axis=1).sort_values(ascending=False)\
        .iloc[0:ntop].index
    df = df.loc[top].transpose()

    x = StandardScaler().fit_transform(df)

    pca = PCA(n_components=n_components)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'], index=df.index)

    if metadata is not None:
        #meta = metadata.loc[principalDf.index]
        #finalDf = pd.concat([principalDf, metadata], axis = 1)
        finalDf = principalDf.merge(metadata, how='left', left_index=True, right_index=True)
    else:
        finalDf = principalDf
    output = dict({
        'df': finalDf,
        'explained_variance_ratio': pca.explained_variance_ratio_
    })
        
    return output

def multiplot_from_generator(g, num_columns, figsize_for_one_row=None):
    # call 'next(g)' to get past the first 'yield'
    next(g)
    # default to 15-inch rows, with square subplots
    if figsize_for_one_row is None:
        figsize_for_one_row = (15, 15/num_columns)
    try:
        while True:
            # call plt.figure once per row
            plt.figure(figsize=figsize_for_one_row)
            for col in range(num_columns):
                ax = plt.subplot(1, num_columns, col+1)
                next(g)
    except StopIteration:
        pass

