import json
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import plotly.express as px
import streamlit as st
import boto3
import awswrangler as wr

boto_session: boto3.Session = boto3.Session(
            region_name="us-west-2",
            profile_name="noble-internal",
            )
s3 = boto_session.client("s3")

component_df = wr.s3.read_excel(
    path="s3://research-formulations/ternary-viscosity/processed/categorized_KFv1_parsed_v2.xlsx",
    boto3_session=boto_session,
)
chem_db = wr.s3.read_csv(
    path="s3://research-formulations/ternary-viscosity/processed/chem_db.csv",
    boto3_session=boto_session,
)
all_form_data = wr.s3.read_csv(
    path="s3://research-formulations/ternary-viscosity/processed/all_form_data.csv",
    boto3_session=boto_session,
)

for col in ['IDs', 'Proportions']:
    all_form_data[col] = all_form_data[col].apply(json.loads)

unique_categories = sorted(set([category for row in component_df[['category_1', 'category_2', 'category_3']].apply(tuple, axis=1) for category in row if type(category) is str]))


st.write("##### Please select a chemical")
if st.toggle("Filter by Category"):
    categories_to_use = st.multiselect(
        "Select Category", 
        unique_categories,
        default=None,
    )
else:
    categories_to_use = unique_categories

def get_filtered_smiles(df, categories):
    return df[df[['category_1', 'category_2', 'category_3']].apply(lambda row: any(category in categories for category in row), axis=1)]['SMILES'].unique().tolist()

smiles_to_use = [""] + get_filtered_smiles(component_df, categories_to_use)

col1, col2 = st.columns((0.7, 0.3))

with col1:
    smi = st.selectbox(
        '',
        smiles_to_use,
        key='select1'
    )

tab1, tab2, tab3 = st.tabs(["Pure Component", "Binary Mixtures", "Ternary or Higher Mixtures"])

if smi != "":
    with col2:
        mol = Chem.MolFromSmiles(smi)

        st.image(
            Draw.MolToImage(mol, size=(600, 600)),
            use_column_width=True,
            caption=f"First component: {smi}",
            width=200,
        )


    with tab1:
        chem_id = chem_db[chem_db['SMILES'] == smi].index[0]

        forms_with_chem = all_form_data[all_form_data['IDs'].apply(lambda x: chem_id in x)]

        # handle formulations with pure component

        forms_with_pure_chem = forms_with_chem[forms_with_chem['Proportions'].apply(lambda x: len(x) == 1)]
        forms_with_pure_chem = forms_with_pure_chem[['T', 'logV']].sort_values(by="T")
        forms_with_pure_chem['T'] = forms_with_pure_chem['T'].round(1)
        forms_with_pure_chem = forms_with_pure_chem.groupby("T").mean().reset_index()
        forms_with_pure_chem.rename(columns={"T": "Temperature (K)", "logV": 'log(Viscosity (cP))'}, inplace=True)
        st.write(f"{len(forms_with_pure_chem)} sample(s) with pure {smi}")
        if len(forms_with_pure_chem) > 1:
            # plotly chart of viscosity vs temperature
            fig = px.line(forms_with_pure_chem, x="Temperature (K)", y="log(Viscosity (cP))")
            fig.update_traces(mode="markers+lines")
            st.plotly_chart(fig)
        else:
            st.dataframe(forms_with_pure_chem, hide_index=True)

    # handle formulations with binary mixtures
    with tab2:
        binary_mixtures_with_chem = forms_with_chem[forms_with_chem['Proportions'].apply(lambda x: len(x) == 2)]
        if len(binary_mixtures_with_chem) == 0:
            st.write(f"No binary mixture(s) containing {smi}")
        if len(binary_mixtures_with_chem) > 0:

            chem_ids_mixed_with = list(set([chem_id_ for forms in binary_mixtures_with_chem['IDs'] for chem_id_ in forms if chem_id_ != chem_id]))
            smis_mixed_with = chem_db.loc[chem_ids_mixed_with]['SMILES'].values

            st.write(f"{len(binary_mixtures_with_chem)} sample(s) with binary mixtures containing {smi}, involving {len(smis_mixed_with)} other chemicals")

            component_df_with_smis_mixed_with = component_df[component_df['SMILES'].isin(smis_mixed_with)]
            categories_mixed_with = sorted(set([category for row in component_df_with_smis_mixed_with[['category_1', 'category_2', 'category_3']].apply(tuple, axis=1) for category in row if type(category) is str]))

            st.write(f'##### Please select a chemical to be mixed with {smi}')
            if st.toggle("Filter by Category", key="filter2"):
                categories_to_use_second_component = st.multiselect(
                    "Select Category", 
                    categories_mixed_with,
                    default=None,
                )
            else:
                categories_to_use_second_component = categories_mixed_with

            smis_mixed_with_filtered = [''] + sorted(set(get_filtered_smiles(component_df, categories_to_use_second_component)) & set(smis_mixed_with))

            col1, col2 = st.columns((0.7, 0.3))

            with col1:
                smi2 = st.selectbox(
                    '',
                    smis_mixed_with_filtered,
                    key='select2',
                )

            if smi2 != "":
                with col2:
                    mol = Chem.MolFromSmiles(smi2)

                    st.image(
                        Draw.MolToImage(mol, size=(600, 600)),
                        use_column_width=True,
                        caption=f"Second component: {smi2}",
                        width=200,
                    )

                chem_id2 = chem_db[chem_db['SMILES'] == smi2].index[0]

                forms_with_both_chem = all_form_data[all_form_data['IDs'].apply(lambda x: chem_id in x and chem_id2 in x)]

                unique_Ts = forms_with_both_chem['T'].unique()
                T = st.selectbox(
                    "Select Temperature",
                    unique_Ts,
                )
                selected_forms_with_both_chems = forms_with_both_chem[forms_with_both_chem['T'] == T]

                def get_amount(row, chem_id):
                    if chem_id in row['IDs']:
                        return row['Proportions'][row['IDs'].index(chem_id)]
                    else:
                        return 0
                    
                selected_forms_with_both_chems[f'{smi} fraction'] = selected_forms_with_both_chems.apply(lambda row: get_amount(row, chem_id), axis=1).round(2)

                selected_forms_with_both_chems = selected_forms_with_both_chems.groupby(f'{smi} fraction')[['logV']].mean().reset_index()
                selected_forms_with_both_chems.rename(columns={"logV": 'log(Viscosity (cP))'}, inplace=True)

                if len(forms_with_pure_chem) > 1:
                    fig = px.line(selected_forms_with_both_chems, x=f'{smi} fraction', y="log(Viscosity (cP))")
                    fig.update_traces(mode="markers+lines")
                    st.plotly_chart(fig)
                else:
                    st.dataframe(selected_forms_with_both_chems[[f'{smi} fraction', "log(Viscosity (cP))"]], hide_index=True)


    with tab3:
        higher_mixtures_with_chem = forms_with_chem[forms_with_chem['Proportions'].apply(lambda x: len(x) > 2)]

        def format_row(row):
            row['IDs'] = ', '.join([chem_db.loc[chem_id]['SMILES'] for chem_id in row['IDs']])
            row['Proportions'] = ', '.join([f"{x:.2g}" for x in row['Proportions']])

            return row
        higher_mixtures_with_chem = higher_mixtures_with_chem.apply(format_row, axis=1).rename(columns={"T": "Temperature (K)", "logV": 'log(Viscosity (cP))'})

        st.write(f"{len(higher_mixtures_with_chem)} sample(s) with ternary or higher mixtures containing {smi}")

        if len(higher_mixtures_with_chem) > 0:
            st.dataframe(higher_mixtures_with_chem, hide_index=True)
