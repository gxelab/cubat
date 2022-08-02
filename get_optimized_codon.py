import statsmodels.api as sm
from genecode_data import *

def optimized_codon(processed_array, enc_array,genecode_data):

    #call genecode database
    
    codon_familyami = genecode_data.codon_familyami
    all_amifamilies = genecode_data.all_amifamilies


    #do the correlation analysis
    
    fop_opt=np.zeros(18,dtype=int)
    zvalue_array=np.zeros(59,dtype=float)
    cir=0
    for i in range(0,len(codon_familyami)-1):
        one_ami= processed_array[:, codon_familyami[i]:codon_familyami[i + 1]]
        ami_count = np.sum(one_ami, axis=1)
        ami_count[ami_count == 0] = 1
        probability=one_ami.T/ami_count.T
        binomial_array=np.zeros(6,dtype=float)
        for m in range(0,all_amifamilies[i]):
            endog= probability[m]
            exog=enc_array
            exog = sm.add_constant(exog, prepend=True)
            res = sm.GLM(endog=endog, var_weights=ami_count, exog=exog,missing=True,
                                    family=sm.families.Binomial()).fit()
            binomial_array[m]=res.tvalues[1]
            zvalue_array[cir]=res.tvalues[1]
            cir+=1
        locate=np.argmin(binomial_array)
        locate+=codon_familyami[i]
        fop_opt[i]=locate
    cbi_opt=np.nonzero(zvalue_array < (-np.sqrt(processed_array.shape[0]) / 3))[0]
    
    return (fop_opt,cbi_opt)
