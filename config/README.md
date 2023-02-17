# Config Files

## Description
These files in a Galaxy installation are used for configuration of the Galaxy instance. Admins can change their specific Galaxy instance to their needs using these.

These files here are put here for reference, the settings and changes made in these files was what made Galaxy able to install the [HybPiper galaxy wrapper](https://github.com/naturalis/galaxy-pipeline-hybseq/tree/main/tools/hybpiper "Hybpiper Galaxy Wrapper folder")'s HybPiper dependency using Conda. 

## galaxy.yml

Usually located at:
*\<Path_To_Your_Galaxy>/config/galaxy.yml*

This file has most of the default options for Galaxy that can be changed here.

The most prominent thing that was changed in this file was the line:

```yml
conda_ensure_channels: iuc,conda-forge,bioconda,chrisjackson-pellicle,defaults
```

## dependency_resolvers_conf.xml

Usually located at:
*\<Path_To_Your_Galaxy>/config/dependency_resolvers_conf.xml* 

This file controls the options for the dependency resolvers in the Galaxy instance. Dependency resolves are where Galaxy installs the requirements for its tools from.

In this file the most prominent option that was changed was:

```xml
<conda />
```
Which was replaced with:
```xml
<conda ensure_channels="conda-forge,bioconda,chrisjackson-pellicle" versionless="false" />
```

The changes in both these files made the chrisjackson-pellicle Conda channel available along with the default Conda Channels.



