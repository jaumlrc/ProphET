
<h4>Installing CPAN:</h4>

Perl modules can be installed using CPAN. Please first certify that CPAN is installed and configured by issuing the command below:

``$ perl -MCPAN -e shell``


* If the command above returns the prompt ```cpan[1]>``` or similar prompt, then CPAN is already configured. So quit the cpan shell by typing: ```cpan[1]> quit```

* If the command returns ```CPAN requires configuration...``` then follow the steps for automatic configuration. Select the default option in every question. Quit CPAN after the configuration is done by typing: ```cpan[1]> quit```

* If the command returns: ```Can't locate CPAN.pm in @INC (@INC contains:... ``` then you will need Administrative privileges to install CPAN either using apt-get ```sudo  apt-get install build-essential``` or yum ```sudo yum install perl-CPAN```

------

<h4>Installing Perl modules:</h4>

Now you are ready to install the required Perl modules. Issue the following commands:

```
$ perl -MCPAN -e 'install Bio::Perl'
$ perl -MCPAN -e 'install LWP::Simple'
$ perl -MCPAN -e 'install XML::Simple'
$ perl -MCPAN -e 'install GD'
```

------
