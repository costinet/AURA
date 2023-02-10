-- phpMyAdmin SQL Dump
-- version 4.0.4
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Generation Time: Aug 11, 2021 at 04:30 PM
-- Server version: 5.6.12-log
-- PHP Version: 5.4.12

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `powerlib`
--

-- --------------------------------------------------------

--
-- Table structure for table `transistor_table`
--

CREATE TABLE IF NOT EXISTS `transistor_table` (
  `partNumber` varchar(50) NOT NULL,
  `type` tinytext,
  `material` tinytext,
  `dateSubmitted` datetime DEFAULT NULL,
  `submittedBy` tinytext,
  `dataSource` varchar(50) DEFAULT NULL,
  `manufacturer` varchar(50) DEFAULT NULL,
  `errorReported` tinyint(1) DEFAULT NULL,
  `width` decimal(3,3) DEFAULT NULL,
  `length` decimal(3,3) DEFAULT NULL,
  `height` decimal(3,3) DEFAULT NULL,
  `weight` decimal(3,3) DEFAULT NULL,
  `package` varchar(20) DEFAULT NULL,
  `Vds_min` float DEFAULT NULL,
  `Vds_typ` float DEFAULT NULL,
  `Vds_max` float DEFAULT NULL,
  `Vds_pulse_min` float DEFAULT NULL,
  `Vds_pulse_typ` float DEFAULT NULL,
  `Vds_pulse_max` float DEFAULT NULL,
  `Vgs_min` float DEFAULT NULL,
  `Vgs_typ` float DEFAULT NULL,
  `Vgs_max` float DEFAULT NULL,
  `Vth_min` float DEFAULT NULL,
  `Vth_typ` float DEFAULT NULL,
  `Vth_max` float DEFAULT NULL,
  `Rds_min` float DEFAULT NULL,
  `Rds_typ` float DEFAULT NULL,
  `Rds_max` float DEFAULT NULL,
  `Id_min` float DEFAULT NULL,
  `Id_typ` float DEFAULT NULL,
  `Id_max` float DEFAULT NULL,
  `Id_pulse_min` float DEFAULT NULL,
  `Id_pulse_typ` float DEFAULT NULL,
  `Id_pulse_max` float DEFAULT NULL,
  `Idss_min` float DEFAULT NULL,
  `Idss_typ` float DEFAULT NULL,
  `Idss_max` float DEFAULT NULL,
  `Igss_min` float DEFAULT NULL,
  `Igss_typ` float DEFAULT NULL,
  `Igss_max` float DEFAULT NULL,
  `Vsd_min` float DEFAULT NULL,
  `Vsd_typ` float DEFAULT NULL,
  `Vsd_max` float DEFAULT NULL,
  `Is_min` float DEFAULT NULL,
  `Is_typ` float DEFAULT NULL,
  `Is_max` float DEFAULT NULL,
  `Pd_min` float DEFAULT NULL,
  `Pd_typ` float DEFAULT NULL,
  `Pd_max` float DEFAULT NULL,
  `Tj_min` float DEFAULT NULL,
  `Tj_typ` float DEFAULT NULL,
  `Tj_max` float DEFAULT NULL,
  `Tstg_min` float DEFAULT NULL,
  `Tstg_typ` float DEFAULT NULL,
  `Tstg_max` float DEFAULT NULL,
  `Tc_min` float DEFAULT NULL,
  `Tc_typ` float DEFAULT NULL,
  `Tc_max` float DEFAULT NULL,
  `Qg_min` float DEFAULT NULL,
  `Qg_typ` float DEFAULT NULL,
  `Qg_max` float DEFAULT NULL,
  `Qgs_min` float DEFAULT NULL,
  `Qgs_typ` float DEFAULT NULL,
  `Qgs_max` float DEFAULT NULL,
  `Qgd_min` float DEFAULT NULL,
  `Qgd_typ` float DEFAULT NULL,
  `Qgd_max` float DEFAULT NULL,
  `Qg_th_min` float DEFAULT NULL,
  `Qg_th_typ` float DEFAULT NULL,
  `Qg_th_max` float DEFAULT NULL,
  `Qoss_min` float DEFAULT NULL,
  `Qoss_typ` float DEFAULT NULL,
  `Qoss_max` float DEFAULT NULL,
  `Qrr_min` float DEFAULT NULL,
  `Qrr_typ` float DEFAULT NULL,
  `Qrr_max` float DEFAULT NULL,
  `Rg_min` float DEFAULT NULL,
  `Rg_typ` float DEFAULT NULL,
  `Rg_max` float DEFAULT NULL,
  `Yfs_min` float DEFAULT NULL,
  `Yfs_typ` float DEFAULT NULL,
  `Yfs_max` float DEFAULT NULL,
  `Eoss_min` float DEFAULT NULL,
  `Eoss_typ` float DEFAULT NULL,
  `Eoss_max` float DEFAULT NULL,
  `Eon_min` float DEFAULT NULL,
  `Eon_typ` float DEFAULT NULL,
  `Eon_max` float DEFAULT NULL,
  `Eoff_min` float DEFAULT NULL,
  `Eoff_typ` float DEFAULT NULL,
  `Eoff_max` float DEFAULT NULL,
  `Ciss_min` float DEFAULT NULL,
  `Ciss_typ` float DEFAULT NULL,
  `Ciss_max` float DEFAULT NULL,
  `Coss_min` float DEFAULT NULL,
  `Coss_typ` float DEFAULT NULL,
  `Coss_max` float DEFAULT NULL,
  `Crss_min` float DEFAULT NULL,
  `Crss_typ` float DEFAULT NULL,
  `Crss_max` float DEFAULT NULL,
  `tr_min` float DEFAULT NULL,
  `tr_typ` float DEFAULT NULL,
  `tr_max` float DEFAULT NULL,
  `tf_min` float DEFAULT NULL,
  `tf_typ` float DEFAULT NULL,
  `tf_max` float DEFAULT NULL,
  `td_off_min` float DEFAULT NULL,
  `td_off_typ` float DEFAULT NULL,
  `td_off_max` float DEFAULT NULL,
  `td_on_min` float DEFAULT NULL,
  `td_on_typ` float DEFAULT NULL,
  `td_on_max` float DEFAULT NULL,
  `Rtheta_jc_min` float DEFAULT NULL,
  `Rtheta_jc_typ` float DEFAULT NULL,
  `Rtheta_jc_max` float DEFAULT NULL,
  `Rtheta_jb_min` float DEFAULT NULL,
  `Rtheta_jb_typ` float DEFAULT NULL,
  `Rtheta_jb_max` float DEFAULT NULL,
  `Rtheta_ja_min` float DEFAULT NULL,
  `Rtheta_ja_typ` float DEFAULT NULL,
  `Rtheta_ja_max` float DEFAULT NULL,
  UNIQUE KEY `partNumber` (`partNumber`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Dumping data for table `transistor_table`
--

INSERT INTO `transistor_table` (`partNumber`, `type`, `material`, `dateSubmitted`, `submittedBy`, `dataSource`, `manufacturer`, `errorReported`, `width`, `length`, `height`, `weight`, `package`, `Vds_min`, `Vds_typ`, `Vds_max`, `Vds_pulse_min`, `Vds_pulse_typ`, `Vds_pulse_max`, `Vgs_min`, `Vgs_typ`, `Vgs_max`, `Vth_min`, `Vth_typ`, `Vth_max`, `Rds_min`, `Rds_typ`, `Rds_max`, `Id_min`, `Id_typ`, `Id_max`, `Id_pulse_min`, `Id_pulse_typ`, `Id_pulse_max`, `Idss_min`, `Idss_typ`, `Idss_max`, `Igss_min`, `Igss_typ`, `Igss_max`, `Vsd_min`, `Vsd_typ`, `Vsd_max`, `Is_min`, `Is_typ`, `Is_max`, `Pd_min`, `Pd_typ`, `Pd_max`, `Tj_min`, `Tj_typ`, `Tj_max`, `Tstg_min`, `Tstg_typ`, `Tstg_max`, `Tc_min`, `Tc_typ`, `Tc_max`, `Qg_min`, `Qg_typ`, `Qg_max`, `Qgs_min`, `Qgs_typ`, `Qgs_max`, `Qgd_min`, `Qgd_typ`, `Qgd_max`, `Qg_th_min`, `Qg_th_typ`, `Qg_th_max`, `Qoss_min`, `Qoss_typ`, `Qoss_max`, `Qrr_min`, `Qrr_typ`, `Qrr_max`, `Rg_min`, `Rg_typ`, `Rg_max`, `Yfs_min`, `Yfs_typ`, `Yfs_max`, `Eoss_min`, `Eoss_typ`, `Eoss_max`, `Eon_min`, `Eon_typ`, `Eon_max`, `Eoff_min`, `Eoff_typ`, `Eoff_max`, `Ciss_min`, `Ciss_typ`, `Ciss_max`, `Coss_min`, `Coss_typ`, `Coss_max`, `Crss_min`, `Crss_typ`, `Crss_max`, `tr_min`, `tr_typ`, `tr_max`, `tf_min`, `tf_typ`, `tf_max`, `td_off_min`, `td_off_typ`, `td_off_max`, `td_on_min`, `td_on_typ`, `td_on_max`, `Rtheta_jc_min`, `Rtheta_jc_typ`, `Rtheta_jc_max`, `Rtheta_jb_min`, `Rtheta_jb_typ`, `Rtheta_jb_max`, `Rtheta_ja_min`, `Rtheta_ja_typ`, `Rtheta_ja_max`) VALUES
('EPC2023', NULL, NULL, '2021-08-11 16:27:50', 'dcostine', 'datasheet', NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 30, NULL, NULL, NULL, -4, NULL, 6, NULL, NULL, NULL, NULL, 1.15, 1.45, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL),
('EPC2040', NULL, NULL, '2021-08-11 16:27:50', 'dcostine', 'datasheet', NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL),
('EPC8002', NULL, NULL, '2021-08-11 16:27:50', 'dcostine', 'datasheet', NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 65, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0.38, 0.48, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
